/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


#include "acados/ocp_qp/ocp_qp_hpmpc.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

acados_size_t ocp_qp_hpmpc_opts_calculate_size(void *config_, ocp_qp_dims *dims)
{
    acados_size_t size = sizeof(ocp_qp_hpmpc_opts);
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



void ocp_qp_hpmpc_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    ocp_qp_hpmpc_opts *opts = opts_;

    if (!strcmp(field, "iter_max"))
    {
        int *tmp_ptr = value;
        opts->max_iter = *tmp_ptr;
    }
    else if (!strcmp(field, "tol_stat"))
    {
        // TODO set solver exit tolerance
    }
    else if (!strcmp(field, "tol_eq"))
    {
        // TODO set solver exit tolerance
    }
    else if (!strcmp(field, "tol_ineq"))
    {
        // TODO set solver exit tolerance
    }
    else if (!strcmp(field, "tol_comp"))
    {
        double *tmp_ptr = value;
        opts->tol = *tmp_ptr;
    }
    else if (!strcmp(field, "warm_start"))
    {
        int *tmp_ptr = value;
        opts->warm_start = *tmp_ptr;
    }
    else
    {
        printf("\nerror: ocp_qp_hpmpc_opts_set: wrong field: %s\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_qp_hpmpc_memory_calculate_size(void *config_, ocp_qp_dims *dims, void *opts_)
{
    ocp_qp_hpmpc_opts *args = (ocp_qp_hpmpc_opts *) opts_;

    size_t N = dims->N;
    size_t N2 = args->N2;
    size_t M = args->M;
    size_t *nx = dims->nx;
    size_t *nu = dims->nu;
    size_t *nb = dims->nb;
    size_t *ng = dims->ng;

    size_t ii;
    acados_size_t ws_size = sizeof(ocp_qp_hpmpc_memory);

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



void ocp_qp_hpmpc_memory_get(void *config_, void *mem_, const char *field, void* value)
{
    qp_solver_config *config = config_;
    ocp_qp_hpmpc_memory *mem = mem_;

    if(!strcmp(field, "time_qp_solver_call"))
    {
        double *tmp_ptr = value;
        *tmp_ptr = mem->time_qp_solver_call;
    }
    else if(!strcmp(field, "iter"))
    {
        int *tmp_ptr = value;
        *tmp_ptr = mem->iter;
    }
    else if (!strcmp(field, "status"))
    {
        int *tmp_ptr = value;
        *tmp_ptr = mem->status;
    }
    else
    {
        printf("\nerror: ocp_qp_hpipm_memory_get: field %s not available\n", field);
        exit(1);
    }

    return;

}

void ocp_qp_hpmpc_memory_reset(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    ocp_qp_in *qp_in = qp_in_;
    // reset memory
    printf("acados: reset hpmpc_mem not implemented.\n");
    exit(1);
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_qp_hpmpc_workspace_calculate_size(void *config_, ocp_qp_dims *dims, void *opts_)
{
    return 0;
}

/************************************************
 * functions
 ************************************************/

int ocp_qp_hpmpc(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    ocp_qp_in *qp_in = qp_in_;
    ocp_qp_out *qp_out = qp_out_;

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
            exit(1);
        }
    }

    qp_info *info = (qp_info *) qp_out->misc;
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
        mem->status = d_ip2_res_mpc_hard_libstr(
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
        mem->status = d_ip2_res_mpc_hard_libstr(
            &kk, k_max, mu0, mu_tol, hpmpc_args->alpha_min, warm_start, mem->stats, N, nx, nu, nb,
            qp_in->idxb, ng, qp_in->BAbt, qp_in->RSQrq, qp_in->DCt, qp_in->d, qp_out->ux,
            compute_mult, mem->hpi, qp_out->lam, qp_out->t,
            mem->hpmpc_work);  // recover original stage M
    }

    info->solve_QP_time = acados_toc(&qp_timer);
    acados_tic(&interface_timer);

    mem->out_iter = kk;  // TODO(dimitris): obsolete

    mem->time_qp_solver_call = info->solve_QP_time;
    mem->iter = kk;

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

    int acados_status = mem->status;
    if (mem->status == 0) acados_status = ACADOS_SUCCESS;
    if (mem->status == 1) acados_status = ACADOS_MAXITER;
    if (mem->status == 2) acados_status = ACADOS_MINSTEP;

    return acados_status;
}


void ocp_qp_hpmpc_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2)
{
    printf("\nerror: ocp_qp_hpmpc_solver_get: not implemented yet\n");
    exit(1);
}

void ocp_qp_hpmpc_eval_sens(void *config_, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: ocp_qp_hpmpc_eval_sens: not implemented yet\n");
    exit(1);
}



void ocp_qp_hpmpc_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->dims_set = &ocp_qp_dims_set;
    config->opts_calculate_size = (size_t (*)(void *, void *)) & ocp_qp_hpmpc_opts_calculate_size;
    config->opts_assign = (void *(*) (void *, void *, void *) ) & ocp_qp_hpmpc_opts_assign;
    config->opts_initialize_default =
        (void (*)(void *, void *, void *)) & ocp_qp_hpmpc_opts_initialize_default;
    config->opts_update = (void (*)(void *, void *, void *)) & ocp_qp_hpmpc_opts_update;
    config->opts_set = &ocp_qp_hpmpc_opts_set;
    config->memory_calculate_size =
        (size_t (*)(void *, void *, void *)) & ocp_qp_hpmpc_memory_calculate_size;
    config->memory_assign =
        (void *(*) (void *, void *, void *, void *) ) & ocp_qp_hpmpc_memory_assign;
    config->memory_get = &ocp_qp_hpmpc_memory_get;
    config->workspace_calculate_size =
        (size_t (*)(void *, void *, void *)) & ocp_qp_hpmpc_workspace_calculate_size;
    config->evaluate = &ocp_qp_hpmpc;
    config->eval_sens = &ocp_qp_hpmpc_eval_sens;
    config->memory_reset = &ocp_qp_hpmpc_memory_reset;
    config->solver_get = &ocp_qp_hpmpc_solver_get;


    return;
}
