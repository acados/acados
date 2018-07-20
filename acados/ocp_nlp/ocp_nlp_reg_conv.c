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

#include <assert.h>
#include <stdlib.h>

#include "acados/utils/math.h"
#include "acados/utils/mem.h"

#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"

int ocp_nlp_reg_conv_calculate_memory_size(ocp_nlp_reg_dims *dims)
{
    int nx = dims->nx, nu = dims->nu;

    int size = 0;

    size += sizeof(ocp_nlp_reg_conv_memory);

    size += nu*nu*sizeof(double);             // R
    size += nu*nu*sizeof(double);             // V
    size += nu*sizeof(double);                // d
    size += (nx+nu)*(nx+nu)*sizeof(double);   // reg_hess

    size += 1 * 64;

    size += blasfeo_memsize_dmat(nx, nx);     // Q_tilde
    size += blasfeo_memsize_dmat(nx, nx);     // Q_bar
    size += blasfeo_memsize_dmat(nx+nu, nx);  // BAQ
    size += blasfeo_memsize_dmat(nu, nu);     // L
    size += blasfeo_memsize_dmat(nx, nx);     // delta_eye
    size += blasfeo_memsize_dmat(nx, nu);     // St_copy

    return size;
}

void *ocp_nlp_reg_conv_assign_memory(ocp_nlp_reg_dims *dims, void *raw_memory)
{
    int nx = dims->nx, nu = dims->nu;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_reg_conv_memory *mem = (ocp_nlp_reg_conv_memory *) c_ptr;

    c_ptr += sizeof(ocp_nlp_reg_conv_memory);

    mem->R = (double *) c_ptr;
    c_ptr += nu*nu*sizeof(double);

    mem->V = (double *) c_ptr;
    c_ptr += nu*nu*sizeof(double);

    mem->d = (double *) c_ptr;
    c_ptr += nu*sizeof(double);

    mem->reg_hess = (double *) c_ptr;
    c_ptr += (nx+nu)*(nx+nu)*sizeof(double);

    align_char_to(64, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->Q_tilde, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->Q_bar, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx+nu, nx, &mem->BAQ, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nu, nu, &mem->L, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->delta_eye, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nu, &mem->St_copy, &c_ptr);

    assert((char *)mem + ocp_nlp_reg_conv_calculate_memory_size(dims) >= c_ptr);

    return mem;
}

void ocp_nlp_reg_conv(void *config, ocp_nlp_reg_dims *dims, ocp_nlp_reg_in *in,
                      ocp_nlp_reg_out *out, ocp_nlp_reg_opts *opts, void *mem_)
{
    int N = dims->N;

    int nx = dims->nx, nu = dims->nu;

    double delta = opts->delta;

    ocp_nlp_reg_conv_memory *mem = (ocp_nlp_reg_conv_memory *) mem_;

    // Algorithm 6 from Verschueren2017

    blasfeo_ddiare(nx, delta, &mem->delta_eye, 0, 0);
    blasfeo_dgecp(nx, nx, &mem->delta_eye, 0, 0, &mem->Q_tilde, 0, 0);
    blasfeo_dgecp(nx, nx, &in->RSQrq[N], 0, 0, &mem->Q_bar, 0, 0);
    blasfeo_dgead(nx, nx, -1.0, &mem->Q_tilde, 0, 0, &mem->Q_bar, 0, 0);

    for (int i = N-1; i >= 0; --i)
    {

        // printf("----------------\n");
        // printf("--- stage %d ---\n", i);
        // printf("----------------\n");

        // printf("QSR\n");
        // blasfeo_print_dmat(nx+nu+1, nx+nu, &work->qp_in->RSQrq[i], 0, 0);

        // printf("Q_bar\n");
        // blasfeo_print_dmat(nx, nx, &Q_bar, 0, 0);

        // printf("BAbt\n");
        // blasfeo_print_dmat(nx+nu, nx, &work->qp_in->BAbt[i], 0, 0);

        blasfeo_dgemm_nt(nx+nu, nx, nx, 1.0, &in->BAbt[i], 0, 0, &mem->Q_bar, 0, 0, 0.0,
                         &mem->BAQ, 0, 0, &mem->BAQ, 0, 0);

        blasfeo_dsyrk_ln_mn(nx+nu+1, nx+nu, nx, 1.0, &in->BAbt[i], 0, 0, &mem->BAQ, 0, 0, 1.0,
                            &in->RSQrq[i], 0, 0, &in->RSQrq[i], 0, 0);

        // printf("BAQ\n");
        // blasfeo_print_dmat(nx+nu, nx, &BAQ, 0, 0);

        blasfeo_unpack_dmat(nu, nu, &in->RSQrq[i], 0, 0, mem->R, nu);
        eigen_decomposition(nu, mem->R, mem->V, mem->d);

        bool needs_regularization = false;
        for (int j = 0; j < nu; ++j)
            if (mem->d[j] < 1e-10)
                needs_regularization = true;

        if (needs_regularization)
        {
            blasfeo_unpack_dmat(nx+nu, nx+nu, &in->RSQrq[i], 0, 0, mem->reg_hess, nx+nu);
            regularize(nx+nu, mem->reg_hess);
            blasfeo_pack_dmat(nx+nu, nx+nu, mem->reg_hess, nx+nu, &in->RSQrq[i], 0, 0);
        }

        // struct blasfeo_dvec *lam = &nlp_out->lam[i];

        // printf("multipliers: %f, %f, diff: %f\n", BLASFEO_DVECEL(lam, 0),
        //        BLASFEO_DVECEL(lam, dims->ni[i]),
        //        fabs(BLASFEO_DVECEL(lam, 0) - BLASFEO_DVECEL(lam, dims->ni[i])));

        // if (fabs(BLASFEO_DVECEL(lam, 0) - BLASFEO_DVECEL(lam, dims->ni[i])) > ACADOS_EPS)
        // {
        //     // active bound

        // }

        // printf("QSR_hat\n");
        // blasfeo_print_dmat(nx+nu+1, nx+nu, &work->qp_in->RSQrq[i], 0, 0);


        blasfeo_dgecp(nx, nx, &in->RSQrq[i], nu, nu, &mem->Q_bar, 0, 0);
        blasfeo_dgecp(nx, nu, &in->RSQrq[i], nu, 0, &mem->St_copy, 0, 0);

        // R = L * L^T
        blasfeo_dpotrf_l(nu, &in->RSQrq[i], 0, 0, &mem->L, 0, 0);
        // Q = S^T * L^-T
        blasfeo_dtrsm_rltn(nx, nu, 1.0, &mem->L, 0, 0, &mem->St_copy, 0, 0, &mem->Q_tilde, 0, 0);

        // Q = S^T * R^-1 * S
        blasfeo_dsyrk_ln(nx, nx, 1.0, &mem->Q_tilde, 0, 0, &mem->Q_tilde, 0, 0,
                         1.0, &mem->delta_eye, 0, 0, &in->RSQrq[i], nu, nu);

        // printf("H_tilde\n");
        // blasfeo_print_dmat(nu+nx, nu+nx, &work->qp_in->RSQrq[i], 0, 0);

        blasfeo_dgead(nx, nx, -1.0, &in->RSQrq[i], nu, nu, &mem->Q_bar, 0, 0);
        blasfeo_dtrtr_l(nx, &mem->Q_bar, 0, 0, &mem->Q_bar, 0, 0);

    }
}
