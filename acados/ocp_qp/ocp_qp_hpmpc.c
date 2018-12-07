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

#include "acados/ocp_qp/ocp_qp_hpmpc.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_v_aux_ext_dep.h"

#include "hpmpc/include/c_interface.h"
#include "hpmpc/include/lqcp_solvers.h"
#include "hpmpc/include/mpc_aux.h"
#include "hpmpc/include/mpc_solvers.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

/************************************************
 * opts
 ************************************************/

int ocp_qp_hpmpc_opts_calculate_size(void *config_, ocp_qp_dims *dims)
{
    int size = sizeof(ocp_qp_hpmpc_opts);
    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;                // align once to typical cache line size
    return size;
}

void *ocp_qp_hpmpc_opts_assign(void *config_, ocp_qp_dims *dims, void *raw_memory)
{
    ocp_qp_hpmpc_opts *args;
    char *c_ptr = (char *) raw_memory;
    args = (ocp_qp_hpmpc_opts *) c_ptr;

    args->N = dims->N;
    args->M = dims->N;
    args->N2 = dims->N;

    c_ptr += sizeof(ocp_qp_hpmpc_opts);

    // align memory to typical cache line size
    size_t s_ptr = (size_t) c_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;
    c_ptr = (char *) s_ptr;

    return (void *) args;
}

void ocp_qp_hpmpc_opts_initialize_default(void *config_, ocp_qp_dims *dims, void *opts_)
{
    ocp_qp_hpmpc_opts *args = (ocp_qp_hpmpc_opts *) opts_;
    args->tol = 1e-8;
    args->max_iter = 100;
    args->mu0 = 1e3;
    args->warm_start = 0;
    args->alpha_min = 1e-8;

    return;
}

void ocp_qp_hpmpc_opts_update(void *config_, ocp_qp_dims *dims, void *opts_)
{
    //    ocp_qp_hpmpc_opts *args = (ocp_qp_hpmpc_opts *)opts_;

    return;
}

/************************************************
 * memory
 ************************************************/

int ocp_qp_hpmpc_memory_calculate_size(void *config_, ocp_qp_dims *dims, void *opts_)
{
    ocp_qp_hpmpc_opts *args = (ocp_qp_hpmpc_opts *) opts_;

    int N = dims->N;
    int N2 = args->N2;
    int M = args->M;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;

    int ii;
    int ws_size = sizeof(ocp_qp_hpmpc_memory);

    ws_size += 6 * args->max_iter * sizeof(double);  // stats

    ws_size += (N + 1) * sizeof(struct blasfeo_dvec);                     // hpi
    for (ii = 0; ii <= N; ii++) ws_size += blasfeo_memsize_dvec(nx[ii]);  // hpi

    ws_size += hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes_noidxb(N, nx, nu, nb, dims->nbx,
                                                                   dims->nbu, ng, N2);

    if (M < N)
    {
        for (ii = 0; ii <= N; ii++)
        {
            ws_size += sizeof(double) * (nu[ii] + nx[ii] + 1) * (nu[ii] + nx[ii]);  // L
            ws_size += sizeof(struct blasfeo_dmat);
            ws_size += sizeof(double) * (nu[ii] + nx[ii]);  // dux
            ws_size += sizeof(struct blasfeo_dvec);
            ws_size += 3 * sizeof(double) * (2 * nb[ii] + 2 * ng[ii]);  // dlam, dt, lamt
            ws_size += 3 * sizeof(struct blasfeo_dvec);
            ws_size += sizeof(double) * (2 * nb[ii] + 2 * ng[ii]);  // tinv
            ws_size += sizeof(struct blasfeo_dvec);
            ws_size += 2 * sizeof(double) * (nb[ii] + ng[ii]);  // Qx, qx
            ws_size += 2 * sizeof(struct blasfeo_dvec);
            ws_size += sizeof(double) * (nx[ii + 1]);  // Pb
            ws_size += sizeof(struct blasfeo_dvec);
        }
        // TODO(all): CHANGE ALL INSTANCES!
        ws_size += blasfeo_memsize_dmat(nx[M] + 1, nx[M]);  // sLxM
        ws_size += sizeof(double) * (nx[M] + 1) * (nx[M]);  // sPpM

        ws_size += d_back_ric_rec_work_space_size_bytes_libstr(N, nx, nu, nb, ng);
    }

    ws_size += 2 * 64;

    return ws_size;
}

void *ocp_qp_hpmpc_memory_assign(void *config_, ocp_qp_dims *dims, void *opts_, void *raw_memory)
{
    ocp_qp_hpmpc_opts *args = (ocp_qp_hpmpc_opts *) opts_;

    // char pointer
    char *c_ptr = (char *) raw_memory;

    ocp_qp_hpmpc_memory *mem = (ocp_qp_hpmpc_memory *) c_ptr;

    c_ptr += sizeof(ocp_qp_hpmpc_memory);

    mem->hpi = (struct blasfeo_dvec *) c_ptr;
    c_ptr += (dims->N + 1) * sizeof(struct blasfeo_dvec);

    mem->stats = (double *) c_ptr;
    c_ptr += 6 * args->max_iter * sizeof(double);

    align_char_to(64, &c_ptr);

    int M = args->M;
    int N = args->N;

    int *nx = (int *) dims->nx;
    int *nu = (int *) dims->nu;
    int *nb = (int *) dims->nb;
    int *ng = (int *) dims->ng;

    int ii;

    if (M < N)
    {
        assign_and_advance_blasfeo_dmat_structs(N + 1, &mem->hsL, &c_ptr);

        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->hsQx, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->hsqx, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->hstinv, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->hsrq, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->hsdux, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->hsdlam, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->hsdpi, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->hsdt, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->hslamt, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->hsPb, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->ux0, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->lam0, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->pi0, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->t0, &c_ptr);

        align_char_to(64, &c_ptr);

        for (ii = 0; ii <= N; ii++)
        {
            // partial tightening-specific
            blasfeo_create_dmat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], &mem->hsL[ii], c_ptr);
            c_ptr += (&mem->hsL[ii])->memsize;
        }

        blasfeo_create_dmat(nx[M] + 1, nx[M], &mem->sLxM, c_ptr);
        c_ptr += (&mem->sLxM)->memsize;
        blasfeo_create_dmat(nx[M] + 1, nx[M], &mem->sPpM, c_ptr);
        c_ptr += (&mem->sPpM)->memsize;

        for (ii = 0; ii <= N; ii++)
        {
            blasfeo_create_dvec(nx[ii] + nu[ii], &mem->ux0[ii], c_ptr);
            c_ptr += (&mem->ux0[ii])->memsize;
            blasfeo_create_dvec(2 * (nb[ii] + ng[ii]), &mem->lam0[ii], c_ptr);
            c_ptr += (&mem->lam0[ii])->memsize;
            blasfeo_create_dvec(2 * (nb[ii] + ng[ii]), &mem->pi0[ii], c_ptr);
            c_ptr += (&mem->pi0[ii])->memsize;
            blasfeo_create_dvec(2 * (nb[ii] + ng[ii]), &mem->t0[ii], c_ptr);
            c_ptr += (&mem->t0[ii])->memsize;

            blasfeo_create_dvec(nu[ii] + nx[ii], &mem->hsdux[ii], c_ptr);
            // blasfeo_pack_dvec(nu[ii]+nx[ii], &mem->ux0[ii], &mem->hsdux[ii], 0);
            c_ptr += (&mem->hsdux[ii])->memsize;

            blasfeo_create_dvec(2 * nb[ii] + 2 * ng[ii], &mem->hstinv[ii], c_ptr);
            c_ptr += (&mem->hstinv[ii])->memsize;
            blasfeo_create_dvec(nb[ii] + ng[ii], &mem->hsQx[ii], c_ptr);
            c_ptr += (&mem->hsQx[ii])->memsize;
            blasfeo_create_dvec(nb[ii] + ng[ii], &mem->hsqx[ii], c_ptr);
            c_ptr += (&mem->hsqx[ii])->memsize;

            blasfeo_create_dvec(2 * nb[ii] + 2 * ng[ii], &mem->hsdlam[ii], c_ptr);
            c_ptr += (&mem->hsdlam[ii])->memsize;
            blasfeo_create_dvec(nx[ii], &mem->hsdpi[ii], c_ptr);
            c_ptr += (&mem->hsdpi[ii])->memsize;
            blasfeo_create_dvec(2 * nb[ii] + 2 * ng[ii], &mem->hsdt[ii], c_ptr);
            c_ptr += (&mem->hsdt[ii])->memsize;
            blasfeo_create_dvec(2 * nb[ii] + 2 * ng[ii], &mem->hslamt[ii], c_ptr);
            c_ptr += (&mem->hslamt[ii])->memsize;

            // partial tightening specific
            blasfeo_create_dvec(nx[ii + 1], &mem->hsPb[ii + 1], c_ptr);
            c_ptr += (&mem->hsPb[ii + 1])->memsize;
        }

        align_char_to(64, &c_ptr);
        mem->work_ric = c_ptr;
        c_ptr += d_back_ric_rec_work_space_size_bytes_libstr(args->N, dims->nx, dims->nu, dims->nb,
                                                             dims->ng);
    }

    for (int ii = 0; ii <= dims->N; ii++)
        assign_and_advance_blasfeo_dvec_mem(dims->nx[ii], &mem->hpi[ii], &c_ptr);

    align_char_to(64, &c_ptr);

    mem->hpmpc_work = (void *) c_ptr;

    // TODO(dimitris): add assert, move hpmpc mem to workspace?
    return raw_memory;
}

/************************************************
 * workspace
 ************************************************/

int ocp_qp_hpmpc_workspace_calculate_size(void *config_, ocp_qp_dims *dims, void *opts_)
{
    return 0;
}

/************************************************
 * functions
 ************************************************/

int ocp_qp_hpmpc(void *config_, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts_, void *mem_,
                 void *work_)
{
    ocp_qp_hpmpc_opts *hpmpc_args = (ocp_qp_hpmpc_opts *) opts_;
    ocp_qp_hpmpc_memory *mem = (ocp_qp_hpmpc_memory *) mem_;

    int N = qp_in->dim->N;
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng;
    int *ns = qp_in->dim->ns;

    for (int ii = 0; ii <= N; ii++)
    {
        if (ns[ii] > 0)
        {
            printf("\nHPMPC interface can not handle ns>0 yet: what about implementing it? :)\n");
            return ACADOS_FAILURE;
        }
    }

    ocp_qp_info *info = (ocp_qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer;
    acados_tic(&tot_timer);

    int ii;

    int hpmpc_status = -1;

    int M = hpmpc_args->M;

    // extract args from struct
    double mu_tol = hpmpc_args->tol;
    int k_max = hpmpc_args->max_iter;
    double mu0 = hpmpc_args->mu0;
    int warm_start = hpmpc_args->warm_start;

    //  other solver arguments
    int compute_mult = 1;
    int kk = -1;  // actual number of iterations

    acados_tic(&interface_timer);

    struct blasfeo_dmat *hsmatdummy = NULL;
    struct blasfeo_dvec *hsvecdummy = NULL;

    for (int ii = 0; ii <= N; ++ii)
    {
        // temporarily invert sign of upper bounds
        blasfeo_dvecsc(nb[ii], -1.0, &qp_in->d[ii], nb[ii] + ng[ii]);
        blasfeo_dvecsc(ng[ii], -1.0, &qp_in->d[ii], 2 * nb[ii] + ng[ii]);
    }

    real_t sigma_mu = hpmpc_args->sigma_mu;

    int nuM;
    int nbM;

    struct blasfeo_dmat hstmpmat0;

    info->interface_time = acados_toc(&interface_timer);
    acados_tic(&qp_timer);

    if (M < N)
    {
        for (int ii = 0; ii <= N; ii++)
            blasfeo_create_dvec(nu[ii] + nx[ii], &mem->hsrq[ii], qp_in->rqz[ii].pa);

        // update cost function matrices and vectors (box constraints)
        d_update_hessian_gradient_mpc_hard_libstr(N - M, &nx[M], &nu[M], &nb[M], &ng[M],
                                                  &qp_in->d[M], sigma_mu, &mem->t0[M],
                                                  &mem->hstinv[M], &mem->lam0[M], &mem->hslamt[M],
                                                  &mem->hsdlam[M], &mem->hsQx[M], &mem->hsqx[M]);

        // backward riccati factorization and solution at the end
        d_back_ric_rec_sv_back_libstr(
            N - M, &nx[M], &nu[M], &nb[M], &qp_in->idxb[M], &ng[M], 0, &qp_in->BAbt[M], hsvecdummy,
            1, &qp_in->RSQrq[M], &mem->hsrq[M], &qp_in->DCt[M], &mem->hsQx[M], &mem->hsqx[M],
            &qp_out->ux[M], 1, &mem->hsdpi[M], 1, &mem->hsPb[M], &mem->hsL[M], mem->work_ric);

        // extract chol factor of [P p; p' *]
        blasfeo_dtrcp_l(nx[M], &mem->hsL[M], nu[M], nu[M], &mem->sLxM, 0, 0);
        blasfeo_dgecp(1, nx[M], &mem->hsL[M], nu[M] + nx[M], nu[M], &mem->sLxM, nx[M], 0);

        // recover [P p; p' *]
        blasfeo_dsyrk_ln_mn(nx[M] + 1, nx[M], nx[M], 1.0, &mem->sLxM, 0, 0, &mem->sLxM, 0, 0, 0.0,
                            &mem->sPpM, 0, 0, &mem->sPpM, 0, 0);

        // backup stage M
        nuM = nu[M];
        nbM = nb[M];
        hstmpmat0 = qp_in->RSQrq[M];

        // update new terminal cost
        nu[M] = 0;
        nb[M] = 0;
        qp_in->RSQrq[M] = mem->sPpM;
        qp_out->ux[M].pa += nuM;

        // IPM at the beginning
        hpmpc_status = d_ip2_res_mpc_hard_libstr(
            &kk, k_max, mu0, mu_tol, hpmpc_args->alpha_min, warm_start, mem->stats, M, nx, nu, nb,
            qp_in->idxb, ng, qp_in->BAbt, qp_in->RSQrq, qp_in->DCt, qp_in->d, qp_out->ux,
            compute_mult, mem->hpi, qp_out->lam, mem->t0,
            mem->hpmpc_work);  // recover original stage M

        nu[M] = nuM;
        nb[M] = nbM;
        qp_in->RSQrq[M] = hstmpmat0;
        qp_out->ux[M].pa -= nuM;

        // forward riccati solution at the end
        d_back_ric_rec_sv_forw_libstr(
            N - M, &nx[M], &nu[M], &nb[M], &qp_in->idxb[M], &ng[M], 0, &qp_in->BAbt[M], hsvecdummy,
            1, &qp_in->RSQrq[M], &mem->hsrq[M], hsmatdummy, &mem->hsQx[M], &mem->hsqx[M],
            &qp_out->ux[M], 1, &mem->hpi[M], 1, &mem->hsPb[M], &mem->hsL[M], mem->work_ric);

        // compute alpha, dlam and dt real_t alpha = 1.0;
        // compute primal step hsdux for stages M to N
        real_t *temp_p1, *temp_p2;
        for (int i = M; i <= N; i++)
        {
            // hsdux is initialized to be equal to hpmpc_args.ux0
            temp_p1 = mem->hsdux[i].pa;
            temp_p2 = qp_out->ux[i].pa;  // hsux[i].pa;
            for (int j = 0; j < nx[i] + nu[i]; j++) temp_p1[j] = -temp_p1[j] + temp_p2[j];
        }

        double alpha;
        d_compute_alpha_mpc_hard_libstr(N - M, &nx[M], &nu[M], &nb[M], &qp_in->idxb[M], &ng[M],
                                        &alpha, &mem->t0[M], &mem->hsdt[M], &qp_out->lam[M],
                                        &mem->hsdlam[M], &mem->hslamt[M], &mem->hsdux[M],
                                        &qp_in->DCt[M], &qp_in->d[M]);

        // update stages M to N
        double mu_scal = 0.0;
        d_update_var_mpc_hard_libstr(N - M, &nx[M], &nu[M], &nb[M], &ng[M], &sigma_mu, mu_scal,
                                     alpha, &qp_out->ux[M], &mem->hsdux[M], &mem->t0[M],
                                     &mem->hsdt[M], &qp_out->lam[M], &mem->hsdlam[M], &mem->hpi[M],
                                     &mem->hsdpi[M]);
    }
    else
    {
        // IPM at the beginning
        hpmpc_status = d_ip2_res_mpc_hard_libstr(
            &kk, k_max, mu0, mu_tol, hpmpc_args->alpha_min, warm_start, mem->stats, N, nx, nu, nb,
            qp_in->idxb, ng, qp_in->BAbt, qp_in->RSQrq, qp_in->DCt, qp_in->d, qp_out->ux,
            compute_mult, mem->hpi, qp_out->lam, qp_out->t,
            mem->hpmpc_work);  // recover original stage M
    }

    info->solve_QP_time = acados_toc(&qp_timer);
    acados_tic(&interface_timer);

    mem->out_iter = kk;  // TODO(dimitris): obsolete

    // copy result to qp_out
    for (ii = 0; ii < N; ii++) blasfeo_dveccp(nx[ii + 1], &mem->hpi[ii + 1], 0, &qp_out->pi[ii], 0);

    // restore sign of upper bounds
    for (int jj = 0; jj <= N; jj++)
    {
        blasfeo_dvecsc(nb[jj], -1.0, &qp_in->d[jj], nb[jj] + ng[jj]);
        blasfeo_dvecsc(ng[jj], -1.0, &qp_in->d[jj], 2 * nb[jj] + ng[jj]);
    }

    info->interface_time += acados_toc(&interface_timer);
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = kk;

    int acados_status = hpmpc_status;
    if (hpmpc_status == 0) acados_status = ACADOS_SUCCESS;
    if (hpmpc_status == 1) acados_status = ACADOS_MAXITER;
    if (hpmpc_status == 2) acados_status = ACADOS_MINSTEP;

    return acados_status;
}

void ocp_qp_hpmpc_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->dims_set = &ocp_qp_dims_set;
    config->opts_calculate_size = (int (*)(void *, void *)) & ocp_qp_hpmpc_opts_calculate_size;
    config->opts_assign = (void *(*) (void *, void *, void *) ) & ocp_qp_hpmpc_opts_assign;
    config->opts_initialize_default =
        (void (*)(void *, void *, void *)) & ocp_qp_hpmpc_opts_initialize_default;
    config->opts_update = (void (*)(void *, void *, void *)) & ocp_qp_hpmpc_opts_update;
    config->memory_calculate_size =
        (int (*)(void *, void *, void *)) & ocp_qp_hpmpc_memory_calculate_size;
    config->memory_assign =
        (void *(*) (void *, void *, void *, void *) ) & ocp_qp_hpmpc_memory_assign;
    config->workspace_calculate_size =
        (int (*)(void *, void *, void *)) & ocp_qp_hpmpc_workspace_calculate_size;
    config->evaluate = (int (*)(void *, void *, void *, void *, void *, void *)) & ocp_qp_hpmpc;

    return;
}
