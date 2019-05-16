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

#include "acados/ocp_nlp/ocp_nlp_reg_convexify.h"

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
 * opts
 ************************************************/

int ocp_nlp_reg_convexify_opts_calculate_size(void)
{
    return sizeof(ocp_nlp_reg_convexify_opts);
}



void *ocp_nlp_reg_convexify_opts_assign(void *raw_memory)
{
    return raw_memory;
}



void ocp_nlp_reg_convexify_opts_initialize_default(void *config_, ocp_nlp_reg_dims *dims, void *opts_)
{
    ocp_nlp_reg_convexify_opts *opts = opts_;

    opts->delta = 1e-4;
    opts->epsilon = 1e-4;

    return;
}



void ocp_nlp_reg_convexify_opts_set(void *config_, ocp_nlp_reg_dims *dims, void *opts_, char *field, void* value)
{

    ocp_nlp_reg_convexify_opts *opts = opts_;

    if (!strcmp(field, "delta"))
    {
        double *d_ptr = value;
        opts->delta = *d_ptr;
    }
    else if (!strcmp(field, "epsilon"))
    {
        double *d_ptr = value;
        opts->epsilon = *d_ptr;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_reg_convexify_opts_set\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * memory
 ************************************************/

int ocp_nlp_reg_convexify_calculate_memory_size(void *config_, ocp_nlp_reg_dims *dims, void *opts_)
{

    int *nx = dims->nx;
    int *nu = dims->nu;
    int N = dims->N;

    int ii;

    int nuM = nu[0];
    for(ii=1; ii<=N; ii++)
    {
        nuM = nu[ii]>nuM ? nu[ii] : nuM;
    }

    int nxM = nx[0];
    for(ii=1; ii<=N; ii++)
    {
        nxM = nx[ii]>nxM ? nx[ii] : nxM;
    }

    int nuxM = nu[0]+nx[0];
    for(ii=1; ii<=N; ii++)
    {
        nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
    }

    int size = 0;

    size += sizeof(ocp_nlp_reg_convexify_memory);

    size += nuM*nuM*sizeof(double);     // R
    size += nuxM*nuxM*sizeof(double);   // V
    size += 2*nuxM*sizeof(double);      // d e
    size += nuxM*nuxM*sizeof(double);   // reg_hess

    size += (N+1)*sizeof(struct blasfeo_dmat); // original_RSQrq

    size += 1 * 64;

    size += 3*blasfeo_memsize_dmat(nxM, nxM);     // Q_tilde Q_bar delta_eye
    size += blasfeo_memsize_dmat(nuxM, nxM);    // BAQ
    size += blasfeo_memsize_dmat(nuM, nuM);     // L
    size += blasfeo_memsize_dmat(nxM, nuM);     // St_copy

    for (ii=0; ii<=N; ii++)
    {
        size += blasfeo_memsize_dmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // original_RSQrq
    }

//    size += 2*blasfeo_memsize_dvec(nxM); // grad b2

    // giaf's
    size += ((N+1)+N)*sizeof(struct blasfeo_dmat *); // RSQrq BAbt
    size += (2*(N+1)+2*N)*sizeof(struct blasfeo_dvec *); // rq b ux pi

    return size;
}



void *ocp_nlp_reg_convexify_assign_memory(void *config_, ocp_nlp_reg_dims *dims, void *opts_, void *raw_memory)
{

    int *nx = dims->nx;
    int *nu = dims->nu;
    int N = dims->N;

    int ii;

    int nuM = nu[0];
    for(ii=1; ii<=N; ii++)
    {
        nuM = nu[ii]>nuM ? nu[ii] : nuM;
    }

    int nxM = nx[0];
    for(ii=1; ii<=N; ii++)
    {
        nxM = nx[ii]>nxM ? nx[ii] : nxM;
    }

    int nuxM = nu[0]+nx[0];
    for(ii=1; ii<=N; ii++)
    {
        nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
    }


    char *c_ptr = (char *) raw_memory;

    ocp_nlp_reg_convexify_memory *mem = (ocp_nlp_reg_convexify_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_reg_convexify_memory);

    mem->R = (double *) c_ptr;
    c_ptr += nuM*nuM*sizeof(double);

    mem->V = (double *) c_ptr;
    c_ptr += nuxM*nuxM*sizeof(double);

    mem->d = (double *) c_ptr;
    c_ptr += nuxM*sizeof(double);

    mem->e = (double *) c_ptr;
    c_ptr += nuxM*sizeof(double);

    mem->reg_hess = (double *) c_ptr;
    c_ptr += nuxM*nuxM*sizeof(double);

    mem->original_RSQrq = (struct blasfeo_dmat *) c_ptr;
    c_ptr += (N+1)*sizeof(struct blasfeo_dmat);

    mem->RSQrq = (struct blasfeo_dmat **) c_ptr;
    c_ptr += (N+1)*sizeof(struct blasfeo_dmat *);

    mem->BAbt = (struct blasfeo_dmat **) c_ptr;
    c_ptr += N*sizeof(struct blasfeo_dmat *);

    mem->rq = (struct blasfeo_dvec **) c_ptr;
    c_ptr += (N+1)*sizeof(struct blasfeo_dvec *);

    mem->b = (struct blasfeo_dvec **) c_ptr;
    c_ptr += N*sizeof(struct blasfeo_dvec *);

    mem->ux = (struct blasfeo_dvec **) c_ptr;
    c_ptr += (N+1)*sizeof(struct blasfeo_dvec *);

    mem->pi = (struct blasfeo_dvec **) c_ptr;
    c_ptr += N*sizeof(struct blasfeo_dvec *);

    align_char_to(64, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nxM, nxM, &mem->Q_tilde, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nxM, nxM, &mem->Q_bar, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nuxM, nxM, &mem->BAQ, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nuM, nuM, &mem->L, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nxM, nxM, &mem->delta_eye, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nxM, nuM, &mem->St_copy, &c_ptr);

    for (ii=0; ii<=N; ii++)
    {
        assign_and_advance_blasfeo_dmat_mem(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &mem->original_RSQrq[ii], &c_ptr);
    }

//    assign_and_advance_blasfeo_dvec_mem(nxM, &mem->grad, &c_ptr);
//    assign_and_advance_blasfeo_dvec_mem(nxM, &mem->b2, &c_ptr);

    assert((char *)mem + ocp_nlp_reg_convexify_calculate_memory_size(config_, dims, opts_) >= c_ptr);

    return mem;
}



void ocp_nlp_reg_convexify_memory_set_RSQrq_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *RSQrq, void *memory_)
{
    ocp_nlp_reg_convexify_memory *memory = memory_;

    int ii;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;

    for(ii=0; ii<=N; ii++)
    {
        memory->RSQrq[ii] = RSQrq+ii;
//        blasfeo_print_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], memory->RSQrq[ii], 0, 0);
    }

    return;
}



void ocp_nlp_reg_convexify_memory_set_rq_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *rq, void *memory_)
{
    ocp_nlp_reg_convexify_memory *memory = memory_;

    int ii;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;

    for(ii=0; ii<=N; ii++)
    {
        memory->rq[ii] = rq+ii;
//        blasfeo_print_dvec(nu[ii]+nx[ii], memory->rq[ii], 0);
    }

    return;
}



void ocp_nlp_reg_convexify_memory_set_BAbt_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *BAbt, void *memory_)
{
    ocp_nlp_reg_convexify_memory *memory = memory_;

    int ii;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;

    for(ii=0; ii<N; ii++)
    {
        memory->BAbt[ii] = BAbt+ii;
//        blasfeo_print_dmat(nu[ii]+nx[ii]+1, nx[ii+1], memory->BAbt[ii], 0, 0);
    }

    return;
}



void ocp_nlp_reg_convexify_memory_set_b_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *b, void *memory_)
{
    ocp_nlp_reg_convexify_memory *memory = memory_;

    int ii;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;

    for(ii=0; ii<N; ii++)
    {
        memory->b[ii] = b+ii;
//        blasfeo_print_dvec(nx[ii=1], memory->b[ii], 0);
    }

    return;
}



void ocp_nlp_reg_convexify_memory_set_ux_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *ux, void *memory_)
{
    ocp_nlp_reg_convexify_memory *memory = memory_;

    int ii;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;

    for(ii=0; ii<=N; ii++)
    {
        memory->ux[ii] = ux+ii;
//        blasfeo_print_dvec(nu[ii]+nx[ii], memory->ux[ii], 0);
    }

    return;
}



void ocp_nlp_reg_convexify_memory_set_pi_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *pi, void *memory_)
{
    ocp_nlp_reg_convexify_memory *memory = memory_;

    int ii;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;

    for(ii=0; ii<N; ii++)
    {
        memory->pi[ii] = pi+ii;
//        blasfeo_print_dvec(nx[ii+1], memory->pi[ii], 0);
    }

    return;
}



void ocp_nlp_reg_convexify_memory_set(void *config_, ocp_nlp_reg_dims *dims, void *memory_, char *field, void *value)
{

    if(!strcmp(field, "RSQrq_ptr"))
    {
        struct blasfeo_dmat *RSQrq = value;
        ocp_nlp_reg_convexify_memory_set_RSQrq_ptr(dims, RSQrq, memory_);
    }
    else if(!strcmp(field, "rq_ptr"))
    {
        struct blasfeo_dvec *rq = value;
        ocp_nlp_reg_convexify_memory_set_rq_ptr(dims, rq, memory_);
    }
    else if(!strcmp(field, "BAbt_ptr"))
    {
        struct blasfeo_dmat *BAbt = value;
        ocp_nlp_reg_convexify_memory_set_BAbt_ptr(dims, BAbt, memory_);
    }
    else if(!strcmp(field, "b_ptr"))
    {
        struct blasfeo_dvec *b = value;
        ocp_nlp_reg_convexify_memory_set_b_ptr(dims, b, memory_);
    }
    else if(!strcmp(field, "ux_ptr"))
    {
        struct blasfeo_dvec *ux = value;
        ocp_nlp_reg_convexify_memory_set_ux_ptr(dims, ux, memory_);
    }
    else if(!strcmp(field, "pi_ptr"))
    {
        struct blasfeo_dvec *pi = value;
        ocp_nlp_reg_convexify_memory_set_pi_ptr(dims, pi, memory_);
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_reg_convexify_set\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * functions
 ************************************************/

// NOTE this only considers the case of (dynamcs) equality constraints (no inequality constraints)
// TODO inequality constraints case
void ocp_nlp_reg_convexify_regularize_hessian(void *config, ocp_nlp_reg_dims *dims, void *opts_, void *mem_)
{
    ocp_nlp_reg_convexify_memory *mem = mem_;
    ocp_nlp_reg_convexify_opts *opts = opts_;

    int ii, jj;

    int *nx = dims->nx;
    int *nu = dims->nu;
    int N = dims->N;

    double delta = opts->delta;

    // Algorithm 6 from Verschueren2017

    blasfeo_dgecp(nu[N]+nx[N]+1, nu[N]+nx[N], mem->RSQrq[N], 0, 0, &mem->original_RSQrq[N], 0, 0);

    // TODO regularize R at last stage if needed !!!
    blasfeo_dgese(nx[N], nu[N]+nx[N], 0.0, &mem->delta_eye, 0, 0);
    blasfeo_ddiare(nx[N], delta, &mem->delta_eye, 0, 0);
    blasfeo_dgecp(nx[N], nx[N], &mem->delta_eye, 0, 0, &mem->Q_tilde, 0, 0);
    blasfeo_dgecp(nx[N], nx[N], mem->RSQrq[N], 0, 0, &mem->Q_bar, 0, 0);
    blasfeo_dgead(nx[N], nx[N], -1.0, &mem->Q_tilde, 0, 0, &mem->Q_bar, 0, 0);
    blasfeo_dtrtr_l(nx[N], &mem->Q_bar, 0, 0, &mem->Q_bar, 0, 0);

    for (ii = N-1; ii >= 0; --ii)
    {
        blasfeo_drowin(nx[ii+1], 1.0, mem->b[ii], 0, mem->BAbt[ii], nu[ii]+nx[ii], 0);
        blasfeo_drowin(nu[ii]+nx[ii], 1.0, mem->rq[ii], 0, mem->RSQrq[ii], nu[ii]+nx[ii], 0);

        blasfeo_dgecp(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0, &mem->original_RSQrq[ii], 0, 0);

        // printf("----------------\n");
        // printf("--- stage %d ---\n", i);
        // printf("----------------\n");

        // printf("QSR\n");
        // blasfeo_print_dmat(nx+nu+1, nx+nu, &work->qp_in->RSQrq[i], 0, 0);

        // printf("Q_bar\n");
        // blasfeo_print_dmat(nx, nx, &Q_bar, 0, 0);

        // printf("BAbt\n");
        // blasfeo_print_dmat(nx+nu, nx, &work->qp_in->BAbt[i], 0, 0);

        // TODO implement using cholesky
        blasfeo_dgemm_nt(nu[ii]+nx[ii], nx[ii], nx[ii+1], 1.0, mem->BAbt[ii], 0, 0, &mem->Q_bar, 0, 0, 0.0, &mem->BAQ, 0, 0, &mem->BAQ, 0, 0);

        blasfeo_dsyrk_ln_mn(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], nx[ii+1], 1.0, mem->BAbt[ii], 0, 0, &mem->BAQ, 0, 0, 1.0, mem->RSQrq[ii], 0, 0, mem->RSQrq[ii], 0, 0);

        // printf("BAQ\n");
        // blasfeo_print_dmat(nx+nu, nx, &BAQ, 0, 0);

        blasfeo_unpack_dmat(nu[ii], nu[ii], mem->RSQrq[ii], 0, 0, mem->R, nu[ii]);
        acados_eigen_decomposition(nu[ii], mem->R, mem->V, mem->d, mem->e);

        bool needs_regularization = false;
        for (jj = 0; jj < nu[ii]; jj++)
            if (mem->d[jj] < 1e-10)
                needs_regularization = true;

        if (needs_regularization)
        {
            blasfeo_unpack_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0, mem->reg_hess, nu[ii]+nx[ii]);
            acados_mirror(nu[ii]+nx[ii], mem->reg_hess, mem->V, mem->d, mem->e, 1e-4);
            blasfeo_pack_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], mem->reg_hess, nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0);
        }

        // printf("QSR_hat\n");
        // blasfeo_print_dmat(nx+nu+1, nx+nu, &work->qp_in->RSQrq[i], 0, 0);


        blasfeo_dgecp(nx[ii], nx[ii], mem->RSQrq[ii], nu[ii], nu[ii], &mem->Q_bar, 0, 0);
        blasfeo_dgecp(nx[ii], nu[ii], mem->RSQrq[ii], nu[ii], 0, &mem->St_copy, 0, 0);

        // R = L * L^T
        blasfeo_dpotrf_l(nu[ii], mem->RSQrq[ii], 0, 0, &mem->L, 0, 0);
        // Q = S^T * L^-T
        blasfeo_dtrsm_rltn(nx[ii], nu[ii], 1.0, &mem->L, 0, 0, &mem->St_copy, 0, 0, &mem->Q_tilde, 0, 0);

        // Q = S^T * R^-1 * S
        blasfeo_dsyrk_ln(nx[ii], nx[ii], 1.0, &mem->Q_tilde, 0, 0, &mem->Q_tilde, 0, 0, 1.0, &mem->delta_eye, 0, 0, mem->RSQrq[ii], nu[ii], nu[ii]);

        // printf("H_tilde\n");
        // blasfeo_print_dmat(nu+nx, nu+nx, &work->qp_in->RSQrq[i], 0, 0);

        // make symmetric
//        blasfeo_dtrtr_l(nx[ii], &mem->Q_bar, 0, 0, &mem->Q_bar, 0, 0);

        // TODO take from b !!!!!!
//        for (jj = 0; jj < nx[ii+1]; jj++)
//            BLASFEO_DVECEL(&mem->b2, jj) = BLASFEO_DMATEL(mem->BAbt[ii], nu[ii]+nx[ii], jj);

        // TODO nx stage is not consistent with above !!!!!!!
//        blasfeo_dgemv_n(nx[ii+1], nx[ii+1], 1.0, &mem->Q_bar, 0, 0, &mem->b2, 0, 0.0, &mem->grad, 0, &mem->grad, 0);
//        blasfeo_dgemv_n(nu[ii]+nx[ii], nx[ii+1], 1.0, mem->BAbt[ii], 0, 0, &mem->grad, 0, 0.0, &mem->b2, 0, &mem->b2, 0);

//        for (jj = 0; jj < nu[ii]+nx[ii]; jj++)
            // TODO maybe 'b' is a bad naming...
//            BLASFEO_DMATEL(mem->RSQrq[ii], nu[ii]+nx[ii], jj) = BLASFEO_DMATEL(mem->RSQrq[ii], nu[ii]+nx[ii], jj) + BLASFEO_DVECEL(&mem->b2, jj);

        blasfeo_dgead(nx[ii], nx[ii], -1.0, mem->RSQrq[ii], nu[ii], nu[ii], &mem->Q_bar, 0, 0);

        // make symmetric
        blasfeo_dtrtr_l(nx[ii], &mem->Q_bar, 0, 0, &mem->Q_bar, 0, 0);

    }

    return;
}



// NOTE this only considers the case of (dynamcs) equality constraints (no inequality constraints)
// TODO inequality constraints case
void ocp_nlp_reg_convexify_correct_dual_sol(void *config, ocp_nlp_reg_dims *dims, void *opts_, void *mem_)
{
    ocp_nlp_reg_convexify_memory *mem = mem_;
    ocp_nlp_reg_convexify_opts *opts = opts_;

    int ii;

    int *nx = dims->nx;
    int *nu = dims->nu;
    int N = dims->N;

    blasfeo_dgemv_n(nx[N], nu[N]+nx[N], 1.0, mem->RSQrq[N], nu[N], 0, mem->ux[N], 0, 1.0, mem->rq[N], nu[N], mem->pi[N-1], 0);

    for(ii=1; ii<N; ii++)
    {
        blasfeo_dgemv_n(nx[N-ii], nu[N-ii]+nx[N-ii], 1.0, mem->RSQrq[N-ii], nu[N-ii], 0, mem->ux[N-ii], 0, 1.0, mem->rq[N-ii], nu[N-ii], mem->pi[N-ii-1], 0);
        blasfeo_dgemv_n(nx[N-ii], nx[N-ii+1], 1.0, mem->BAbt[N-ii], nu[N-ii], 0, mem->pi[N-ii], 0, 1.0, mem->pi[N-ii-1], 0, mem->pi[N-ii-1], 0);
    }

    return;
}



void ocp_nlp_reg_convexify_config_initialize_default(ocp_nlp_reg_config *config)
{
    // dims
    config->dims_calculate_size = &ocp_nlp_reg_dims_calculate_size;
    config->dims_assign = &ocp_nlp_reg_dims_assign;
    config->dims_set = &ocp_nlp_reg_dims_set;
    // opts
    config->opts_calculate_size = &ocp_nlp_reg_convexify_opts_calculate_size;
    config->opts_assign = &ocp_nlp_reg_convexify_opts_assign;
    config->opts_initialize_default = &ocp_nlp_reg_convexify_opts_initialize_default;
    config->opts_set = &ocp_nlp_reg_convexify_opts_set;
    // memory
    config->memory_calculate_size = &ocp_nlp_reg_convexify_calculate_memory_size;
    config->memory_assign = &ocp_nlp_reg_convexify_assign_memory;
    config->memory_set = &ocp_nlp_reg_convexify_memory_set;
    config->memory_set_RSQrq_ptr = &ocp_nlp_reg_convexify_memory_set_RSQrq_ptr;
    config->memory_set_rq_ptr = &ocp_nlp_reg_convexify_memory_set_rq_ptr;
    config->memory_set_BAbt_ptr = &ocp_nlp_reg_convexify_memory_set_BAbt_ptr;
    config->memory_set_b_ptr = &ocp_nlp_reg_convexify_memory_set_b_ptr;
    config->memory_set_ux_ptr = &ocp_nlp_reg_convexify_memory_set_ux_ptr;
    config->memory_set_pi_ptr = &ocp_nlp_reg_convexify_memory_set_pi_ptr;
    // functions
    config->regularize_hessian = &ocp_nlp_reg_convexify_regularize_hessian;
    config->correct_dual_sol = &ocp_nlp_reg_convexify_correct_dual_sol;
}
