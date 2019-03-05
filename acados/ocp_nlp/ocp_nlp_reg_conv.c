/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "acados/ocp_nlp/ocp_nlp_reg_conv.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "acados/utils/math.h"
#include "acados/utils/mem.h"

#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"



/************************************************
 * memory
 ************************************************/

int ocp_nlp_reg_conv_calculate_memory_size(void *config_, ocp_nlp_reg_dims *dims, void *opts_)
{
    int nx = dims->nx[0], nu = dims->nu[0], N = dims->N;

    int size = 0;

    size += sizeof(ocp_nlp_reg_conv_memory);

    size += nu*nu*sizeof(double);             // R
    size += (nu+nx)*(nu+nx)*sizeof(double);   // V
    size += 2*(nu+nx)*sizeof(double);         // d e
    size += (nx+nu)*(nx+nu)*sizeof(double);   // reg_hess
    size += (N+1)*sizeof(struct blasfeo_dmat); // original_RSQrq

    size += 1 * 64;

    size += blasfeo_memsize_dmat(nx, nx);     // Q_tilde
    size += blasfeo_memsize_dmat(nx, nx);     // Q_bar
    size += blasfeo_memsize_dmat(nx+nu, nx);  // BAQ
    size += blasfeo_memsize_dmat(nu, nu);     // L
    size += blasfeo_memsize_dmat(nx, nx);     // delta_eye
    size += blasfeo_memsize_dmat(nx, nu);     // St_copy

    size += blasfeo_memsize_dmat(nx+1, nx);
    for (int i = N-1; i >= 0; --i)
    {
        size += blasfeo_memsize_dmat(nu+nx+1, nu+nx);
    }

    size += blasfeo_memsize_dvec(nx);
    size += blasfeo_memsize_dvec(nx);

	size += (N+1)*sizeof(struct blasfeo_dmat *); // RSQrq

    return size;
}



void *ocp_nlp_reg_conv_assign_memory(void *config_, ocp_nlp_reg_dims *dims, void *opts_, void *raw_memory)
{
    int nx = dims->nx[0], nu = dims->nu[0], N = dims->N;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_reg_conv_memory *mem = (ocp_nlp_reg_conv_memory *) c_ptr;

    c_ptr += sizeof(ocp_nlp_reg_conv_memory);

    mem->R = (double *) c_ptr;
    c_ptr += nu*nu*sizeof(double);

    mem->V = (double *) c_ptr;
    c_ptr += (nu+nx)*(nu+nx)*sizeof(double);

    mem->d = (double *) c_ptr;
    c_ptr += (nu+nx)*sizeof(double);

    mem->e = (double *) c_ptr;
    c_ptr += (nu+nx)*sizeof(double);

    mem->reg_hess = (double *) c_ptr;
    c_ptr += (nx+nu)*(nx+nu)*sizeof(double);

    mem->original_RSQrq = (struct blasfeo_dmat *) c_ptr;
    c_ptr += (N+1)*sizeof(struct blasfeo_dmat);

	mem->RSQrq = (struct blasfeo_dmat **) c_ptr;
	c_ptr += (N+1)*sizeof(struct blasfeo_dmat *); // RSQrq

    align_char_to(64, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->Q_tilde, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->Q_bar, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx+nu, nx, &mem->BAQ, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nu, nu, &mem->L, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->delta_eye, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nu, &mem->St_copy, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nx+1, nx, &mem->original_RSQrq[N], &c_ptr);
    for (int i = N-1; i >= 0; --i)
    {
        assign_and_advance_blasfeo_dmat_mem(nu+nx+1, nu+nx, &mem->original_RSQrq[i], &c_ptr);
    }

    assign_and_advance_blasfeo_dvec_mem(nx, &mem->grad, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, &mem->b, &c_ptr);

    assert((char *)mem + ocp_nlp_reg_conv_calculate_memory_size(config_, dims, opts_) >= c_ptr);

    return mem;
}



void ocp_nlp_reg_conv_memory_set_RSQrq_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *RSQrq, void *memory_)
{
    ocp_nlp_reg_conv_memory *memory = memory_;

	int ii;

	int N = dims->N;
	int *nx = dims->nx;
	int *nu = dims->nu;

	for(ii=0; ii<=N; ii++)
	{
		memory->RSQrq[ii] = RSQrq+ii;
//		blasfeo_print_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], memory->RSQrq[ii], 0, 0);
	}

    return;
}



void ocp_nlp_reg_conv_memory_set_BAbt_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *BAbt, void *memory_)
{
    ocp_nlp_reg_conv_memory *memory = memory_;

	int ii;

	int N = dims->N;
	int *nx = dims->nx;
	int *nu = dims->nu;

	for(ii=0; ii<N; ii++)
	{
		memory->BAbt[ii] = BAbt+ii;
//		blasfeo_print_dmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], memory->BAbt[ii], 0, 0);
	}

    return;
}



void ocp_nlp_reg_conv_memory_set(void *config_, ocp_nlp_reg_dims *dims, void *memory_, char *field, void *value)
{

	if(!strcmp(field, "RSQrq_ptr"))
	{
		struct blasfeo_dmat *RSQrq = value;
		ocp_nlp_reg_conv_memory_set_RSQrq_ptr(dims, RSQrq, memory_);
	}
	if(!strcmp(field, "BAbt_ptr"))
	{
		struct blasfeo_dmat *BAbt = value;
		ocp_nlp_reg_conv_memory_set_BAbt_ptr(dims, BAbt, memory_);
	}
	else
	{
		printf("\nerror: field %s not available in ocp_nlp_reg_conv_set\n", field);
		exit(1);
	}

    return;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_reg_conv(void *config, ocp_nlp_reg_dims *dims, void *opts_, void *mem_)
{
    ocp_nlp_reg_conv_memory *mem = mem_;
    ocp_nlp_reg_opts *opts = opts_;

    int N = dims->N;

    int nx = dims->nx[0], nu = dims->nu[0];

    double delta = opts->delta;

    // Algorithm 6 from Verschueren2017

    blasfeo_dgecp(nx+1, nx, mem->RSQrq[N], 0, 0, &mem->original_RSQrq[N], 0, 0);

    blasfeo_ddiare(nx, delta, &mem->delta_eye, 0, 0);
    blasfeo_dgecp(nx, nx, &mem->delta_eye, 0, 0, &mem->Q_tilde, 0, 0);
    blasfeo_dgecp(nx, nx, mem->RSQrq[N], 0, 0, &mem->Q_bar, 0, 0);
    blasfeo_dgead(nx, nx, -1.0, &mem->Q_tilde, 0, 0, &mem->Q_bar, 0, 0);

    for (int i = N-1; i >= 0; --i)
    {
        blasfeo_dgecp(nu+nx+1, nu+nx, mem->RSQrq[i], 0, 0, &mem->original_RSQrq[i], 0, 0);

        // printf("----------------\n");
        // printf("--- stage %d ---\n", i);
        // printf("----------------\n");

        // printf("QSR\n");
        // blasfeo_print_dmat(nx+nu+1, nx+nu, &work->qp_in->RSQrq[i], 0, 0);

        // printf("Q_bar\n");
        // blasfeo_print_dmat(nx, nx, &Q_bar, 0, 0);

        // printf("BAbt\n");
        // blasfeo_print_dmat(nx+nu, nx, &work->qp_in->BAbt[i], 0, 0);

        blasfeo_dgemm_nt(nx+nu, nx, nx, 1.0, mem->BAbt[i], 0, 0, &mem->Q_bar, 0, 0, 0.0,
                         &mem->BAQ, 0, 0, &mem->BAQ, 0, 0);

        blasfeo_dsyrk_ln_mn(nx+nu+1, nx+nu, nx, 1.0, mem->BAbt[i], 0, 0, &mem->BAQ, 0, 0, 1.0,
                            mem->RSQrq[i], 0, 0, mem->RSQrq[i], 0, 0);

        // printf("BAQ\n");
        // blasfeo_print_dmat(nx+nu, nx, &BAQ, 0, 0);

        blasfeo_unpack_dmat(nu, nu, mem->RSQrq[i], 0, 0, mem->R, nu);
        acados_eigen_decomposition(nu, mem->R, mem->V, mem->d, mem->e);

        bool needs_regularization = false;
        for (int j = 0; j < nu; ++j)
            if (mem->d[j] < 1e-10)
                needs_regularization = true;

        if (needs_regularization)
        {
            blasfeo_unpack_dmat(nx+nu, nx+nu, mem->RSQrq[i], 0, 0, mem->reg_hess, nx+nu);
            acados_mirror(nx+nu, mem->reg_hess, mem->V, mem->d, mem->e, 1e-4);
            blasfeo_pack_dmat(nx+nu, nx+nu, mem->reg_hess, nx+nu, mem->RSQrq[i], 0, 0);
        }

        // printf("QSR_hat\n");
        // blasfeo_print_dmat(nx+nu+1, nx+nu, &work->qp_in->RSQrq[i], 0, 0);


        blasfeo_dgecp(nx, nx, mem->RSQrq[i], nu, nu, &mem->Q_bar, 0, 0);
        blasfeo_dgecp(nx, nu, mem->RSQrq[i], nu, 0, &mem->St_copy, 0, 0);

        // R = L * L^T
        blasfeo_dpotrf_l(nu, mem->RSQrq[i], 0, 0, &mem->L, 0, 0);
        // Q = S^T * L^-T
        blasfeo_dtrsm_rltn(nx, nu, 1.0, &mem->L, 0, 0, &mem->St_copy, 0, 0, &mem->Q_tilde, 0, 0);

        // Q = S^T * R^-1 * S
        blasfeo_dsyrk_ln(nx, nx, 1.0, &mem->Q_tilde, 0, 0, &mem->Q_tilde, 0, 0,
                         1.0, &mem->delta_eye, 0, 0, mem->RSQrq[i], nu, nu);

        // printf("H_tilde\n");
        // blasfeo_print_dmat(nu+nx, nu+nx, &work->qp_in->RSQrq[i], 0, 0);

        // make symmetric
        blasfeo_dtrtr_l(nx, &mem->Q_bar, 0, 0, &mem->Q_bar, 0, 0);

        for (int k = 0; k < nx; ++k)
            BLASFEO_DVECEL(&mem->b, k) = BLASFEO_DMATEL(mem->BAbt[i], nu+nx, k);

        blasfeo_dgemv_n(nx, nx, 1.0, &mem->Q_bar, 0, 0, &mem->b, 0, 0.0,
                        &mem->grad, 0, &mem->grad, 0);
        blasfeo_dgemv_n(nu+nx, nx, 1.0, mem->BAbt[i], 0, 0, &mem->grad, 0, 0.0,
                        &mem->b, 0, &mem->b, 0);

        for (int k = 0; k < nu+nx; ++k)
            BLASFEO_DMATEL(mem->RSQrq[i], nu+nx, k) = BLASFEO_DMATEL(mem->RSQrq[i], nu+nx, k)
                                                      + BLASFEO_DVECEL(&mem->b, k);

        blasfeo_dgead(nx, nx, -1.0, mem->RSQrq[i], nu, nu, &mem->Q_bar, 0, 0);
        blasfeo_dtrtr_l(nx, &mem->Q_bar, 0, 0, &mem->Q_bar, 0, 0);

    }
}



void ocp_nlp_reg_conv_config_initialize_default(ocp_nlp_reg_config *config)
{
	// dims
    config->dims_calculate_size = &ocp_nlp_reg_dims_calculate_size;
    config->dims_assign = &ocp_nlp_reg_dims_assign;
    config->dims_set = &ocp_nlp_reg_dims_set;
	// opts
    config->opts_calculate_size = &ocp_nlp_reg_opts_calculate_size;
    config->opts_assign = &ocp_nlp_reg_opts_assign;
    config->opts_initialize_default = &ocp_nlp_reg_opts_initialize_default;
    config->opts_set = &ocp_nlp_reg_opts_set;
	// memory
    config->memory_calculate_size = &ocp_nlp_reg_conv_calculate_memory_size;
    config->memory_assign = &ocp_nlp_reg_conv_assign_memory;
    config->memory_set = &ocp_nlp_reg_conv_memory_set;
    config->memory_set_RSQrq_ptr = &ocp_nlp_reg_conv_memory_set_RSQrq_ptr;
    config->memory_set_BAbt_ptr = &ocp_nlp_reg_conv_memory_set_BAbt_ptr;
	// functions
    config->evaluate = &ocp_nlp_reg_conv;
}
