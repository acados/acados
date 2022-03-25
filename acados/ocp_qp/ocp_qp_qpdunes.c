/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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


// TODO(dimitris): revive linear MPC example

#include "acados/ocp_qp/ocp_qp_qpdunes.h"

// external
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// acados
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"



static int get_maximum_number_of_inequality_constraints(ocp_qp_dims *dims)
{
    int nDmax = dims->ng[0];

    for (int kk = 1; kk < dims->N + 1; kk++)
    {
        if (dims->ng[kk] > nDmax) nDmax = dims->ng[kk];
    }
    return nDmax;
}



static qpdunes_stage_qp_solver_t check_stage_qp_solver(ocp_qp_qpdunes_opts *opts, ocp_qp_in *qp_in)
{
    int N = qp_in->dim->N;
    int nx = qp_in->dim->nx[0];
    int nu = qp_in->dim->nu[0];
    int nD = 0;

    qpdunes_stage_qp_solver_t stageQpSolver = QPDUNES_WITH_CLIPPING;

    // check for polyhedral constraints
    for (int kk = 0; kk < N + 1; kk++) nD += qp_in->dim->ng[kk];
    if (nD != 0) stageQpSolver = QPDUNES_WITH_QPOASES;

    // check for non-zero cross terms
    for (int kk = 0; kk < N; kk++)
    {
        for (int ii = 0; ii < nx; ii++)
        {
            for (int jj = 0; jj < nu; jj++)
            {
                if (BLASFEO_DMATEL(&qp_in->RSQrq[kk], ii + nu, jj) != 0)
                {
                    stageQpSolver = QPDUNES_WITH_QPOASES;
                }
            }
        }
    }

    // check for non-diagonal Q
    int nu_k = nu;
    for (int kk = 0; kk < N + 1; kk++)
    {
        if (kk == N) nu_k = 0;
        for (int ii = 0; ii < nx; ii++)
        {
            for (int jj = 0; jj < nx; jj++)
            {
                if ((ii != jj) && (BLASFEO_DMATEL(&qp_in->RSQrq[kk], ii + nu_k, jj + nu_k) != 0))
                {
                    stageQpSolver = QPDUNES_WITH_QPOASES;
                }
            }
        }
    }

    // check for non-diagonal R
    for (int kk = 0; kk < N; kk++)
    {
        for (int ii = 0; ii < nu; ii++)
        {
            for (int jj = 0; jj < nu; jj++)
            {
                if ((ii != jj) && (BLASFEO_DMATEL(&qp_in->RSQrq[kk], ii, jj) != 0))
                {
                    stageQpSolver = QPDUNES_WITH_QPOASES;
                }
            }
        }
    }
    return stageQpSolver;
}



/************************************************
 * opts
 ************************************************/

acados_size_t ocp_qp_qpdunes_opts_calculate_size(void *config_, ocp_qp_dims *dims)
{
    acados_size_t size = 0;
    size += sizeof(ocp_qp_qpdunes_opts);
    return size;
}



void *ocp_qp_qpdunes_opts_assign(void *config_, ocp_qp_dims *dims, void *raw_memory)
{
    ocp_qp_qpdunes_opts *opts;

    char *c_ptr = (char *) raw_memory;

    opts = (ocp_qp_qpdunes_opts *) c_ptr;
    c_ptr += sizeof(ocp_qp_qpdunes_opts);

    assert((char *) raw_memory + ocp_qp_qpdunes_opts_calculate_size(config_, dims) >= c_ptr);

    return (void *) opts;
}



void ocp_qp_qpdunes_opts_initialize_default(void *config_, ocp_qp_dims *dims, void *opts_)
{
    ocp_qp_qpdunes_opts *opts = (ocp_qp_qpdunes_opts *) opts_;

    // TODO(dimitris): this should be type for all QP solvers and be passed in init. default opts
    qpdunes_options_t qpdunes_opts = QPDUNES_ACADO_SETTINGS;

    opts->stageQpSolver = QPDUNES_WITH_QPOASES;

    opts->options = qpDUNES_setupDefaultOptions();
    opts->isLinearMPC = 0;
    opts->options.printLevel = 0;
    opts->options.stationarityTolerance = 1e-12;
    opts->warmstart = 1;

    if (qpdunes_opts == QPDUNES_DEFAULT_ARGUMENTS)
    {
        // keep default options
    }
    else if (qpdunes_opts == QPDUNES_NONLINEAR_MPC)
    {
        // not implemented yet
    }
    else if (qpdunes_opts == QPDUNES_LINEAR_MPC)
    {
        opts->isLinearMPC = 1;
    }
    else if (qpdunes_opts == QPDUNES_ACADO_SETTINGS)
    {
        opts->options.maxIter = 1000;
        opts->options.printLevel = 0;
        opts->options.stationarityTolerance = 1.e-6;
        opts->options.regParam = 1.e-6;
        opts->options.newtonHessDiagRegTolerance = 1.e-8;
        if (opts->stageQpSolver == QPDUNES_WITH_QPOASES)
            opts->options.lsType = QPDUNES_LS_HOMOTOPY_GRID_SEARCH;
        else if (opts->stageQpSolver == QPDUNES_WITH_CLIPPING)
            opts->options.lsType = QPDUNES_LS_ACCELERATED_GRADIENT_BISECTION_LS;
        opts->options.lineSearchReductionFactor = 0.1;
        opts->options.lineSearchMaxStepSize = 1.;
        opts->options.maxNumLineSearchIterations = 25;
        opts->options.maxNumLineSearchRefinementIterations = 25;
        opts->options.regType = QPDUNES_REG_LEVENBERG_MARQUARDT;
    }
    else
    {
        printf("\nUnknown option (%d) for qpDUNES!\n", qpdunes_opts);
    }

    return;
}



void ocp_qp_qpdunes_opts_update(void *config_, ocp_qp_dims *dims, void *opts_)
{
    //    ocp_qp_qpdunes_opts *opts = (ocp_qp_qpdunes_opts *)opts_;

    return;
}



void ocp_qp_qpdunes_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    ocp_qp_qpdunes_opts *opts = opts_;

    if (!strcmp(field, "tol_stat"))
    {
        double *tol = value;
        // opts->options.stationarityTolerance = *tol;
        // NOTE: there seems to be a mismatch between stationarity tolarance in qpDUNES and acados.
        // SQP didnt converge without the factor 1e-1 here.
        opts->options.stationarityTolerance = 1e-1 * *tol;
    }
    else if (!strcmp(field, "tol_eq"))
    {
        double *tol = value;
        opts->options.equalityTolerance = *tol;
    }
    else if (!strcmp(field, "tol_ineq"))
    {
        double *tol = value;
        opts->options.activenessTolerance = *tol;
    }
    else if (!strcmp(field, "tol_comp"))
    {
        // NOTE(oj): I think qpDUNES does not have this, not sure about though.
    }
    else if (!strcmp(field, "warm_start"))
    {
        // TODO set solver warm start
    }
    else if (!strcmp(field, "iter_max"))
    {
        int *iter_max = value;
        opts->options.maxIter = *iter_max;
    }
    else
    {
        printf("\nerror: ocp_qp_qpdunes_opts_set: wrong field: %s\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_qp_qpdunes_memory_calculate_size(void *config_, ocp_qp_dims *dims, void *opts_)
{
    // NOTE(dimitris): calculate size does NOT include the memory required by qpDUNES
    acados_size_t size = 0;
    size += sizeof(ocp_qp_qpdunes_memory);
    return size;
}



void *ocp_qp_qpdunes_memory_assign(void *config_, ocp_qp_dims *dims, void *opts_, void *raw_memory)
{
    ocp_qp_qpdunes_opts *opts = (ocp_qp_qpdunes_opts *) opts_;
    ocp_qp_qpdunes_memory *mem;

    // char pointer
    char *c_ptr = (char *) raw_memory;

    mem = (ocp_qp_qpdunes_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_qpdunes_memory);

    // initialize memory
    int N, nx, nu;
    unsigned int *nD_ptr = 0;

    N = dims->N;
    nx = dims->nx[0];
    nu = dims->nu[0];

    mem->firstRun = 1;
    mem->nx = nx;
    mem->nu = nu;
    mem->nz = nx + nu;
    mem->nDmax = get_maximum_number_of_inequality_constraints(dims);

    // Check on dimensions
    for (int kk = 1; kk < N; kk++)
    {
        if (dims->nx[kk] != nx || dims->nu[kk] != nu)
        {
            printf("\nocp_qp_qpdunes_memory_assign: nx and nu must be constant for intermediate stages!\n");
            printf("\tGot nx[0] = %d, nu[0] = %d, at stage %d: nx = %d, nu = %d!\n", nx, nu, kk, dims->nx[kk], dims->nu[kk]);
            exit(1);
        }
    }
    if (dims->nx[N] != nx)
    {
        printf("\nocp_qp_qpdunes_memory_assign: nx must be constant for intermediate + terminal stages!\n");
        exit(1);
    }
    if (dims->nu[N] != 0)
    {
        printf("\nocp_qp_qpdunes_memory_assign: nu must be zero for terminal stage!\n");
        exit(1);
    }

    if (opts->stageQpSolver == QPDUNES_WITH_QPOASES)
    {
        // NOTE(dimitris): lsType 5 seems to work but yields WRONG results with ineq. constraints!
        if (opts->options.lsType != 7)
        {
            opts->options.lsType = 7;
            // TODO(dimitris): write proper acados warnings and errors
            if (opts->options.printLevel > 0)
            {
                printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                printf("WARNING: Changed line-search algorithm for qpDUNES (incompatible with QP)");
                printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            }
        }
        if (mem->nDmax > 0)
            nD_ptr = (unsigned int *) dims->ng;  // otherwise leave pointer equal to zero
    }

    // qpDUNES memory allocation
    return_t return_value = qpDUNES_setup(&(mem->qpData), N, nx, nu, nD_ptr, &(opts->options));
    if (return_value != QPDUNES_OK) return NULL;

    return mem;
}



void ocp_qp_qpdunes_memory_get(void *config_, void *mem_, const char *field, void* value)
{
    // qp_solver_config *config = config_;
    ocp_qp_qpdunes_memory *mem = mem_;

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
        printf("\nerror: ocp_qp_qpdunes_memory_get: field %s not available\n", field);
        exit(1);
    }

    return;

}


void ocp_qp_qpdunes_memory_reset(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    // ocp_qp_in *qp_in = qp_in_;
    // reset memory
    printf("acados: reset qpdunes_mem not implemented.\n");
    exit(1);
}



static void form_H(double *H, int nx, int nu, struct blasfeo_dmat *sRSQrq)
{
    // make Q full
    blasfeo_dtrtr_l(nx, sRSQrq, nu, nu, sRSQrq, nu, nu);
    // copy Q
    blasfeo_unpack_dmat(nx, nx, sRSQrq, nu, nu, &H[0], nx + nu);

    // make R full
    blasfeo_dtrtr_l(nu, sRSQrq, 0, 0, sRSQrq, 0, 0);
    // copy R
    blasfeo_unpack_dmat(nu, nu, sRSQrq, 0, 0, &H[nx * (nx + nu) + nx], nx + nu);

    // copy S
    blasfeo_unpack_dmat(nx, nu, sRSQrq, nu, 0, &H[nx * (nx + nu)], nx + nu);
    blasfeo_unpack_tran_dmat(nx, nu, sRSQrq, nu, 0, &H[nx], nx + nu);

    // printf("acados RSQ (nx = %d, nu = %d)\n", nx, nu);
    // blasfeo_print_dmat(sRSQrq->m-1, sRSQrq->n, sRSQrq, 0, 0);
    // printf("qpDUNES Hessian:\n");
    // d_print_mat(nx+nu, nx+nu, H, nx+nu);
    // printf("********************************************\n\n");
}



static void form_RSQ(double *R, double *S, double *Q, int nx, int nu, struct blasfeo_dmat *sRSQrq)
{
    // make Q full
    blasfeo_dtrtr_l(nx, sRSQrq, nu, nu, sRSQrq, nu, nu);
    // copy Q
    blasfeo_unpack_dmat(nx, nx, sRSQrq, nu, nu, Q, nx);

    // make R full
    blasfeo_dtrtr_l(nu, sRSQrq, 0, 0, sRSQrq, 0, 0);
    // copy R
    blasfeo_unpack_dmat(nu, nu, sRSQrq, 0, 0, R, nu);

    // copy S
    blasfeo_unpack_tran_dmat(nx, nu, sRSQrq, nu, 0, S, nu);

    // printf("acados RSQ (nx = %d, nu = %d)\n", nx, nu);
    // blasfeo_print_dmat(sRSQrq->m-1, sRSQrq->n, sRSQrq, 0, 0);
    // printf("qpDUNES Q':\n");
    // d_print_mat(nx, nx, Q, nx);
    // printf("qpDUNES R':\n");
    // d_print_mat(nu, nu, R, nu);
    // printf("qpDUNES S':\n");
    // d_print_mat(nu, nx, S, nu);
    // // d_print_mat(nx, nu, S, nx);
    // printf("********************************************\n\n");
}



static void form_g(double *g, int nx, int nu, struct blasfeo_dvec *srq)
{
    blasfeo_unpack_dvec(nx, srq, nu, &g[0], 1);
    blasfeo_unpack_dvec(nu, srq, 0, &g[nx], 1);

    // printf("acados rq (nx = %d, nu = %d)\n", nx, nu);
    // blasfeo_print_tran_dvec(srq->m, srq, 0);
    // printf("qpDUNES g':\n");
    // d_print_mat(1, nx+nu, g, 1);
    // printf("********************************************\n\n");
}



static void form_dynamics(double *ABt, double *b, int nx, int nu, struct blasfeo_dmat *sBAbt,
                          struct blasfeo_dvec *sb)
{
    // copy A
    blasfeo_unpack_dmat(nx, nx, sBAbt, nu, 0, &ABt[0], nx + nu);
    // copy B
    blasfeo_unpack_dmat(nu, nx, sBAbt, 0, 0, &ABt[nx], nx + nu);
    // copy b
    blasfeo_unpack_dvec(nx, sb, 0, b, 1);

    // printf("acados [B'; A'] (nx = %d, nu = %d)\n", nx, nu);
    // blasfeo_print_dmat(sBAbt->m-1, sBAbt->n, sBAbt, 0, 0);
    // printf("qpDUNES A:\n");
    // d_print_mat(nx, nx, &ABt[0], nx+nu);
    // printf("qpDUNES B:\n");
    // d_print_mat(nu, nx, &ABt[nx], nx+nu);
    // printf("********************************************\n\n");
}



static void form_bounds(double *zLow, double *zUpp, int nx, int nu, int nb, int ng, int *idxb,
                        struct blasfeo_dvec *sd, double infty)
{
    for (int ii = 0; ii < nx + nu; ii++)
    {
        zLow[ii] = -infty;
        zUpp[ii] = infty;
    }
    for (int ii = 0; ii < nb; ii++)
    {
        if (idxb[ii] < nu)
        {                                                             // input bound
            zLow[idxb[ii] + nx] = BLASFEO_DVECEL(sd, ii);             // lb[ii]
            zUpp[idxb[ii] + nx] = -BLASFEO_DVECEL(sd, ii + nb + ng);  // ub[ii]
        }
        else
        {                                                             // state bounds
            zLow[idxb[ii] - nu] = BLASFEO_DVECEL(sd, ii);             // lb[ii]
            zUpp[idxb[ii] - nu] = -BLASFEO_DVECEL(sd, ii + nb + ng);  // ub[ii]
        }
    }
}



static void form_inequalities(double *Ct, double *lc, double *uc, int nx, int nu, int nb, int ng,
                              struct blasfeo_dmat *sDCt, struct blasfeo_dvec *sd)
{
    int ii;

    // copy C
    blasfeo_unpack_dmat(nx, ng, sDCt, nu, 0, &Ct[0], nx + nu);
    // copy D
    blasfeo_unpack_dmat(nu, ng, sDCt, 0, 0, &Ct[nx], nx + nu);
    // copy lc
    blasfeo_unpack_dvec(ng, sd, nb, lc, 1);
    // copy uc
    blasfeo_unpack_dvec(ng, sd, 2 * nb + ng, uc, 1);

    for (ii = 0; ii < ng; ii++) uc[ii] = -uc[ii];
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_qp_qpdunes_workspace_calculate_size(void *config_, ocp_qp_dims *dims, void *opts_)
{
    int nx = dims->nx[0];
    int nu = dims->nu[0];
    int nDmax = get_maximum_number_of_inequality_constraints(dims);
    size_t nz = nx + nu;

    acados_size_t size = sizeof(ocp_qp_qpdunes_workspace);

    size += nz * nz * sizeof(double);     // H
    size += nx * nx * sizeof(double);     // Q
    size += nu * nu * sizeof(double);     // R
    size += nx * nu * sizeof(double);     // S
    size += nz * sizeof(double);          // g
    size += nx * nz * sizeof(double);     // ABt
    size += nz * sizeof(double);          // b
    size += nDmax * nz * sizeof(double);  // Ct
    size += 2 * nDmax * sizeof(double);   // lc, uc
    size += 2 * nz * sizeof(double);      // zLow, zUpp

    return size;
}



static void ocp_qp_qpdunes_cast_workspace(ocp_qp_qpdunes_workspace *work, ocp_qp_qpdunes_memory *mem)
{
    char *c_ptr = (char *) work;
    c_ptr += sizeof(ocp_qp_qpdunes_workspace);

    int nx = mem->nx;
    int nu = mem->nu;
    int nz = mem->nz;
    int nDmax = mem->nDmax;

    assign_and_advance_double(nz * nz, &work->H, &c_ptr);
    assign_and_advance_double(nx * nx, &work->Q, &c_ptr);
    assign_and_advance_double(nu * nu, &work->R, &c_ptr);
    assign_and_advance_double(nx * nu, &work->S, &c_ptr);
    assign_and_advance_double(nz, &work->g, &c_ptr);
    assign_and_advance_double(nx * nz, &work->ABt, &c_ptr);
    assign_and_advance_double(nz, &work->b, &c_ptr);
    assign_and_advance_double(nDmax * nz, &work->Ct, &c_ptr);
    assign_and_advance_double(nDmax, &work->lc, &c_ptr);
    assign_and_advance_double(nDmax, &work->uc, &c_ptr);
    assign_and_advance_double(nz, &work->zLow, &c_ptr);
    assign_and_advance_double(nz, &work->zUpp, &c_ptr);
}



static int update_memory(ocp_qp_in *in, ocp_qp_qpdunes_opts *opts, ocp_qp_qpdunes_memory *mem,
                         ocp_qp_qpdunes_workspace *work)
{
    boolean_t isLTI;  // TODO(dimitris): use isLTI flag for LTI systems
    return_t value = 0;
    qpdunes_stage_qp_solver_t stageQps;

    int N = in->dim->N;
    int nx = in->dim->nx[0];
    int nu = in->dim->nu[0];
    int *nb = in->dim->nb;
    int *ng = in->dim->ng;

    // coldstart
    if (opts->warmstart == 0)
    {
        for (int ii = 0; ii < N; ii++)
            for (int jj = 0; jj < nx; jj++) mem->qpData.lambda.data[ii * nx + jj] = 0.0;
    }

    mem->qpData.options.maxIter = opts->options.maxIter;

    if (mem->firstRun == 1)
    {
        // check if qpDUNES will detect clipping or qpOASES
        stageQps = check_stage_qp_solver(opts, in);

        if (opts->stageQpSolver == QPDUNES_WITH_CLIPPING && stageQps == QPDUNES_WITH_QPOASES)
            return QPDUNES_ERR_INVALID_ARGUMENT;  // user specified clipping but problem requires
                                                  // qpOASES

        // if user specified qpOASES but clipping is detected, trick qpDUNES to detect qpOASES
        // NOTE(dimitris): also needed when partial condensing is used because qpDUNES detects
        // clipping on last stage and crashes..
        if (opts->stageQpSolver == QPDUNES_WITH_QPOASES)
        {
            // make Q[N] non diagonal
            BLASFEO_DMATEL(&in->RSQrq[N], 0, 1) += 1e-8;
            BLASFEO_DMATEL(&in->RSQrq[N], 1, 0) += 1e-8;
        }

        // setup of intervals
        for (int kk = 0; kk < N; ++kk)
        {
            form_RSQ(work->R, work->S, work->Q, nx, nu, &in->RSQrq[kk]);
            form_g(work->g, nx, nu, &in->rqz[kk]);
            form_bounds(work->zLow, work->zUpp, nx, nu, nb[kk], ng[kk], in->idxb[kk], &in->d[kk],
                        opts->options.QPDUNES_INFTY);
            form_dynamics(work->ABt, work->b, nx, nu, &in->BAbt[kk], &in->b[kk]);

            if (opts->stageQpSolver == QPDUNES_WITH_QPOASES)
            {
                if (ng[kk] == 0)
                {
                    value = qpDUNES_setupRegularInterval(&(mem->qpData), mem->qpData.intervals[kk],
                                                         0, work->Q, work->R, work->S, work->g,
                                                         work->ABt, 0, 0, work->b, work->zLow,
                                                         work->zUpp, 0, 0, 0, 0, 0, 0, 0);
                }
                else
                {
                    form_inequalities(work->Ct, work->lc, work->uc, nx, nu, nb[kk], ng[kk],
                                      &in->DCt[kk], &in->d[kk]);
                    value = qpDUNES_setupRegularInterval(
                        &(mem->qpData), mem->qpData.intervals[kk], 0, work->Q, work->R, work->S,
                        work->g, work->ABt, 0, 0, work->b, work->zLow, work->zUpp, 0, 0, 0, 0,
                        work->Ct, work->lc, work->uc);
                }
            }
            else
            {  // do not pass S[kk] or Cx[kk]/Cu[kk] at all
                value = qpDUNES_setupRegularInterval(
                    &(mem->qpData), mem->qpData.intervals[kk], 0, work->Q, work->R, 0, work->g,
                    work->ABt, 0, 0, work->b, work->zLow, work->zUpp, 0, 0, 0, 0, 0, 0, 0);
            }
            if (value != QPDUNES_OK)
            {
                printf("Setup of qpDUNES failed on interval %d\n", kk);
                return (int) value;
            }
        }
        form_bounds(work->zLow, work->zUpp, nx, 0, nb[N], ng[N], in->idxb[N], &in->d[N],
                    opts->options.QPDUNES_INFTY);
        form_RSQ(work->R, work->S, work->Q, nx, 0, &in->RSQrq[N]);
        form_g(work->g, nx, 0, &in->rqz[N]);  // work->g = q
        if (ng[N] == 0)
        {
            value = qpDUNES_setupFinalInterval(&(mem->qpData), mem->qpData.intervals[N], work->Q,
                                               work->g, work->zLow, work->zUpp, 0, 0, 0);
        }
        else
        {
            form_inequalities(work->Ct, work->lc, work->uc, nx, 0, nb[N], ng[N], &in->DCt[N],
                              &in->d[N]);
            value = qpDUNES_setupFinalInterval(&(mem->qpData), mem->qpData.intervals[N], work->Q,
                                               work->g, work->zLow, work->zUpp, work->Ct, work->lc,
                                               work->uc);
        }
        if (value != QPDUNES_OK)
        {
            printf("Setup of qpDUNES failed on last interval\n");
            return (int) value;
        }

        // setup of stage QPs
        value = qpDUNES_setupAllLocalQPs(&(mem->qpData), isLTI = QPDUNES_FALSE);
        if (value != QPDUNES_OK)
        {
            printf("Setup of qpDUNES failed on initialization of stage QPs\n");
            return (int) value;
        }
    }
    else
    {  // if mem->firstRun == 0
        if (opts->isLinearMPC == 0)
        {
            for (int kk = 0; kk < N; kk++)
            {
                form_H(work->H, nx, nu, &in->RSQrq[kk]);
                form_g(work->g, nx, nu, &in->rqz[kk]);
                form_dynamics(work->ABt, work->b, nx, nu, &in->BAbt[kk], &in->b[kk]);

                form_bounds(work->zLow, work->zUpp, nx, nu, nb[kk], ng[kk], in->idxb[kk],
                            &in->d[kk], opts->options.QPDUNES_INFTY);
                if (ng[kk] == 0)
                {
                    value = qpDUNES_updateIntervalData(&(mem->qpData), mem->qpData.intervals[kk],
                                                       work->H, work->g, work->ABt, work->b,
                                                       work->zLow, work->zUpp, 0, 0, 0, 0);
                }
                else
                {
                    form_inequalities(work->Ct, work->lc, work->uc, nx, nu, nb[kk], ng[kk],
                                      &in->DCt[kk], &in->d[kk]);
                    value = qpDUNES_updateIntervalData(
                        &(mem->qpData), mem->qpData.intervals[kk], work->H, work->g, work->ABt,
                        work->b, work->zLow, work->zUpp, work->Ct, work->lc, work->uc, 0);
                }
                if (value != QPDUNES_OK)
                {
                    printf("Update of qpDUNES failed on interval %d\n", kk);
                    return (int) value;
                }
                // qpDUNES_printMatrixData( work->ABt, nx, nx+nu, "AB[%d]", kk);
            }
            form_bounds(work->zLow, work->zUpp, nx, 0, nb[N], ng[N], in->idxb[N], &in->d[N],
                        opts->options.QPDUNES_INFTY);
            form_RSQ(work->R, work->S, work->Q, nx, 0, &in->RSQrq[N]);
            form_g(work->g, nx, 0, &in->rqz[N]);  // work->g = q
            if (ng[N] == 0)
            {
                value =
                    qpDUNES_updateIntervalData(&(mem->qpData), mem->qpData.intervals[N], work->Q,
                                               work->g, 0, 0, work->zLow, work->zUpp, 0, 0, 0, 0);
            }
            else
            {
                form_inequalities(work->Ct, work->lc, work->uc, nx, 0, nb[N], ng[N], &in->DCt[N],
                                  &in->d[N]);
                value = qpDUNES_updateIntervalData(&(mem->qpData), mem->qpData.intervals[N],
                                                   work->Q, work->g, 0, 0, work->zLow, work->zUpp,
                                                   work->Ct, work->lc, work->uc, 0);
            }
            if (value != QPDUNES_OK)
            {
                printf("Update of qpDUNES failed on last interval\n");
                return (int) value;
            }
        }
        else
        {  // linear MPC
            form_bounds(work->zLow, work->zUpp, nx, nu, nb[0], ng[0], in->idxb[0], &in->d[0],
                        opts->options.QPDUNES_INFTY);
            value = qpDUNES_updateIntervalData(&(mem->qpData), mem->qpData.intervals[0], 0, 0, 0, 0,
                                               work->zLow, work->zUpp, 0, 0, 0, 0);
            if (value != QPDUNES_OK)
            {
                printf("Update of qpDUNES failed on first interval\n");
                return (int) value;
            }
        }
    }
    mem->firstRun = 0;
    return (int) value;
}



static void fill_in_qp_out(ocp_qp_in *in, ocp_qp_out *out, ocp_qp_qpdunes_memory *mem)
{
    int N = in->dim->N;
    int *nx = in->dim->nx;
    int *nu = in->dim->nu;
    int *ng = in->dim->ng;
    int *nb = in->dim->nb;

    int **idxb = in->idxb;

    for (int kk = 0; kk < N + 1; kk++)
    {
        int nv = mem->qpData.intervals[kk]->nV;
        int nc = mem->qpData.intervals[kk]->nD;
        double *dual_sol = &mem->qpData.intervals[kk]->y.data[0];
        for (int ii = 0; ii < 2 * nv + 2 * nc; ii++)
            dual_sol[ii] = (dual_sol[ii] >= 0.0) ? dual_sol[ii] : 0.0;

        blasfeo_pack_dvec(nx[kk], &mem->qpData.intervals[kk]->z.data[0], 1, &out->ux[kk], nu[kk]);
        blasfeo_pack_dvec(nu[kk], &mem->qpData.intervals[kk]->z.data[nx[kk]], 1, &out->ux[kk], 0);

        for (int ii = 0; ii < 2 * nb[kk] + 2 * ng[kk]; ii++) out->lam[kk].pa[ii] = 0.0;

        // qpdunes multipliers form pairs (lam_low, lam_upp), first for bounds, then for constraints
        for (int ii = 0; ii < nb[kk]; ii++)
        {
            // flip bounds from [x u] to [u x]
            if (ii < nx[kk])
            {
                out->lam[kk].pa[nu[kk] + ii] = dual_sol[2 * idxb[kk][ii]];
                out->lam[kk].pa[nb[kk] + ng[kk] + nu[kk] + ii] = dual_sol[2 * idxb[kk][ii] + 1];
            }
            else
            {
                out->lam[kk].pa[ii - nx[kk]] = dual_sol[2 * idxb[kk][ii]];
                out->lam[kk].pa[nb[kk] + ng[kk] + ii - nx[kk]] = dual_sol[2 * idxb[kk][ii] + 1];
            }
        }

        for (int ii = 0; ii < ng[kk]; ii++)
        {
            out->lam[kk].pa[nb[kk] + ii] = dual_sol[2 * nv + 2 * ii];
            out->lam[kk].pa[2 * nb[kk] + ng[kk] + ii] = dual_sol[2 * nv + 2 * ii + 1];
        }
    }
    for (int kk = 0; kk < N; kk++)
    {
        blasfeo_pack_dvec(nx[kk + 1], &mem->qpData.lambda.data[kk * nx[kk + 1]], 1, &out->pi[kk], 0);
    }
}



void ocp_qp_qpdunes_free_memory(void *mem_)
{
    ocp_qp_qpdunes_memory *mem = (ocp_qp_qpdunes_memory *) mem_;
    qpDUNES_cleanup(&(mem->qpData));
}



/************************************************
 * functions
 ************************************************/

int ocp_qp_qpdunes(void *config_, ocp_qp_in *in, ocp_qp_out *out, void *opts_, void *mem_,
                   void *work_)
{
    int N = in->dim->N;
    int *ns = in->dim->ns;

    for (int ii = 0; ii <= N; ii++)
    {
        if (ns[ii] > 0)
        {
            printf("\nqpDUNES interface can not handle ns>0 yet: what about implementing it? :)\n");
            exit(1);
        }
    }

    acados_timer tot_timer, qp_timer, interface_timer;
    qp_info *info = (qp_info *) out->misc;
    acados_tic(&tot_timer);

    // cast data structures
    ocp_qp_qpdunes_opts *opts = (ocp_qp_qpdunes_opts *) opts_;
    ocp_qp_qpdunes_memory *mem = (ocp_qp_qpdunes_memory *) mem_;
    ocp_qp_qpdunes_workspace *work = (ocp_qp_qpdunes_workspace *) work_;

    acados_tic(&interface_timer);
    ocp_qp_qpdunes_cast_workspace(work, mem);
    return_t qpdunes_status = update_memory(in, opts, mem, work);
    if (qpdunes_status != QPDUNES_OK) return qpdunes_status;
    info->interface_time = acados_toc(&interface_timer);

    acados_tic(&qp_timer);
    qpdunes_status = qpDUNES_solve(&(mem->qpData));
    info->solve_QP_time = acados_toc(&qp_timer);

    mem->time_qp_solver_call = info->solve_QP_time;
    mem->iter = mem->qpData.log.numIter;

    acados_tic(&interface_timer);
    fill_in_qp_out(in, out, mem);
    ocp_qp_compute_t(in, out);
    info->interface_time += acados_toc(&interface_timer);

    info->total_time = acados_toc(&tot_timer);
    info->num_iter = mem->qpData.log.numIter;
    info->t_computed = 1;

    int acados_status = qpdunes_status;
    if ( qpdunes_status == QPDUNES_OK ||
         qpdunes_status == QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND )
    {
        acados_status = ACADOS_SUCCESS;
    }
    else if (qpdunes_status == QPDUNES_SUCC_SUBOPTIMAL_TERMINATION)
    {
        printf("\nqpDUNES: Success, terminated with suboptimal soluation.\n");
    }
    else if (qpdunes_status == QPDUNES_ERR_ITERATION_LIMIT_REACHED)
    {
        printf("\nqpDUNES: reached maximum number of iterations.\n");
        acados_status = ACADOS_MAXITER;
    }
    else if (qpdunes_status == QPDUNES_ERR_DIVISION_BY_ZERO)
    {
        printf("\nqpDUNES: error, division by zero.\n");
        // NOTE(oj): qpDUNES seems quite conservative here.
        printf("NOTE: the reason might be a not well conditioned hessian.\n");
        printf("qpDUNES seems to have really strong requirements with this regard...\n");
        acados_status = ACADOS_QP_FAILURE;
    }
    else
    {
        printf("\nqpDUNES: returned error not handled by acados.\n");
        acados_status = ACADOS_QP_FAILURE;
    }

    mem->status = acados_status;
    return acados_status;
}



void ocp_qp_qpdunes_eval_sens(void *config_, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: ocp_qp_qpdunes_eval_sens: not implemented yet\n");
    exit(1);
}



void ocp_qp_qpdunes_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->dims_set = &ocp_qp_dims_set;
    config->opts_calculate_size = (acados_size_t (*)(void *, void *)) & ocp_qp_qpdunes_opts_calculate_size;
    config->opts_assign = (void *(*) (void *, void *, void *) ) & ocp_qp_qpdunes_opts_assign;
    config->opts_initialize_default =
        (void (*)(void *, void *, void *)) & ocp_qp_qpdunes_opts_initialize_default;
    config->opts_update = (void (*)(void *, void *, void *)) & ocp_qp_qpdunes_opts_update;
    config->opts_set = &ocp_qp_qpdunes_opts_set;
    config->memory_calculate_size =
        (acados_size_t (*)(void *, void *, void *)) & ocp_qp_qpdunes_memory_calculate_size;
    config->memory_assign =
        (void *(*) (void *, void *, void *, void *) ) & ocp_qp_qpdunes_memory_assign;
    config->memory_get = &ocp_qp_qpdunes_memory_get;
    config->memory_reset = &ocp_qp_qpdunes_memory_reset;
    config->workspace_calculate_size =
        (acados_size_t (*)(void *, void *, void *)) & ocp_qp_qpdunes_workspace_calculate_size;
    config->evaluate = (int (*)(void *, void *, void *, void *, void *, void *)) & ocp_qp_qpdunes;
    config->eval_sens = &ocp_qp_qpdunes_eval_sens;

    return;
}
