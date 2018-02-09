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

// TODO(dimitris): revive linear MPC example

#include "acados/ocp_qp/ocp_qp_qpdunes.h"

//external
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
// acados
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/mem.h"


static int get_maximum_number_of_inequality_constraints(ocp_qp_dims *dims)
{
    int nDmax = dims->ng[0];

    for (int kk = 1; kk < dims->N + 1; kk++)
    {
        if (dims->ng[kk] > nDmax) nDmax = dims->ng[kk];
    }
    return nDmax;
}



static qpdunes_stage_qp_solver_t check_stage_qp_solver(ocp_qp_qpdunes_args *args, ocp_qp_in *qp_in)
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
                if (DMATEL_LIBSTR(&qp_in->RSQrq[kk], ii+nu, jj) != 0)
                {
                    stageQpSolver = QPDUNES_WITH_QPOASES;
                }
            }
        }
    }

    // check for non-diagonal Q
    int nu_k = nu;
    for (int kk = 0; kk < N+1; kk++)
    {
        if (kk == N) nu_k = 0;
        for (int ii = 0; ii < nx; ii++)
        {
            for (int jj = 0; jj < nx; jj++)
            {
                if ((ii != jj) && (DMATEL_LIBSTR(&qp_in->RSQrq[kk], ii+nu_k, jj+nu_k) != 0))
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
                if ((ii != jj) && (DMATEL_LIBSTR(&qp_in->RSQrq[kk], ii, jj) != 0))
                {
                    stageQpSolver = QPDUNES_WITH_QPOASES;
                }
            }
        }
    }
    return stageQpSolver;
}



int ocp_qp_qpdunes_calculate_args_size(ocp_qp_dims *dims, void *submodules_)
{
    int size = 0;
    size += sizeof(ocp_qp_qpdunes_args);
    return size;
}



void *ocp_qp_qpdunes_assign_args(ocp_qp_dims *dims, void **submodules_, void *raw_memory)
{
    ocp_qp_qpdunes_args *args;

    char *c_ptr = (char *) raw_memory;

    args = (ocp_qp_qpdunes_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_qpdunes_args);

    assert((char*)raw_memory + ocp_qp_qpdunes_calculate_args_size(dims, *submodules_) == c_ptr);

    // Update submodules pointer
    *submodules_ = NULL;

    return (void *)args;
}



void *ocp_qp_qpdunes_copy_args(ocp_qp_dims *dims, void *raw_memory, void *source_)
{
    ocp_qp_qpdunes_args *source = (ocp_qp_qpdunes_args *) source_;
    ocp_qp_qpdunes_args *dest;

    void *submodules;

    dest = ocp_qp_qpdunes_assign_args(dims, &submodules, raw_memory);

    dest->options = source->options;
    dest->stageQpSolver = source->stageQpSolver;
    dest->warmstart = source->warmstart;
    dest->isLinearMPC = source->isLinearMPC;

    return (void *)dest;
}



void ocp_qp_qpdunes_initialize_default_args(void *args_)
{
    ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args *)args_;

    // TODO(dimitris): this should be type for all QP solvers and be passed in init. default args
    qpdunes_options_t opts = QPDUNES_NONLINEAR_MPC;

    args->stageQpSolver = QPDUNES_WITH_CLIPPING;

    args->options = qpDUNES_setupDefaultOptions();
    args->isLinearMPC = 0;
    args->options.printLevel = 0;
    args->options.stationarityTolerance = 1e-12;
    args->warmstart = 1;

    if (opts == QPDUNES_DEFAULT_ARGUMENTS) {
    } else if (opts == QPDUNES_NONLINEAR_MPC) {
    } else if (opts == QPDUNES_LINEAR_MPC) {
        args->isLinearMPC = 1;
    } else {
        printf("\nUnknown option (%d) for qpDUNES!\n", opts);
    }
}



int ocp_qp_qpdunes_calculate_memory_size(ocp_qp_dims *dims, void *args_)
{
    // NOTE(dimitris): calculate size does NOT include the memory required by qpDUNES
    int size = 0;
    size += sizeof(ocp_qp_qpdunes_memory);
    return size;
}



void *ocp_qp_qpdunes_assign_memory(ocp_qp_dims *dims, void *args_, void *raw_memory)
{
    ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args *)args_;
    ocp_qp_qpdunes_memory *mem;

    // char pointer
    char *c_ptr = (char *)raw_memory;

    mem = (ocp_qp_qpdunes_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_qpdunes_memory);

    // initialize memory
    int N, nx, nu;
    uint *nD_ptr = 0;
    return_t return_value;

    N = dims->N;
    nx = dims->nx[0];
    nu = dims->nu[0];

    mem->firstRun = 1;
    mem->nx = nx;
    mem->nu = nu;
    mem->nz = nx+nu;
    mem->nDmax = get_maximum_number_of_inequality_constraints(dims);

    // Check for constant dimensions
    for (int kk = 1; kk < N; kk++)
    {
        assert(dims->nx[kk] == nx && "qpDUNES does not support varying dimensions");
        assert(dims->nu[kk] == nu && "qpDUNES does not support varying dimensions");
    }
    assert(dims->nx[N] == nx && "qpDUNES does not support varying dimensions");
    assert(dims->nu[N] == 0 && "qpDUNES does not support nu > 0 on terminal stage");

    if (args->stageQpSolver == QPDUNES_WITH_QPOASES)
    {
        // NOTE(dimitris): lsType 5 seems to work but yields WRONG results with ineq. constraints!
        if (args->options.lsType != 7)
        {
            args->options.lsType = 7;
            // TODO(dimitris): write proper acados warnings and errors
            if (args->options.printLevel > 0)
            {
                printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                printf("WARNING: Changed line-search algorithm for qpDUNES (incompatible with QP)");
                printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            }
        }
        if (mem->nDmax > 0) nD_ptr = (uint *)dims->ng;  // otherwise leave pointer equal to zero
    }

    // qpDUNES memory allocation
    return_value = qpDUNES_setup(&(mem->qpData), N, nx, nu, nD_ptr, &(args->options));
    assert(return_value == QPDUNES_OK && "Setup of the QP solver failed\n");

    return mem;
}



static void form_H(double *H, int nx, int nu, struct blasfeo_dmat *sRSQrq)
{
    // make Q full
    blasfeo_dtrtr_l(nx, sRSQrq, nu, nu, sRSQrq, nu, nu);
    // copy Q
    blasfeo_unpack_dmat(nx, nx, sRSQrq, nu, nu, &H[0], nx+nu);

    // make R full
    blasfeo_dtrtr_l(nu, sRSQrq, 0, 0, sRSQrq, 0, 0);
    // copy R
    blasfeo_unpack_dmat(nu, nu, sRSQrq, 0, 0, &H[nx*(nx+nu) + nx], nx+nu);

    // copy S
    blasfeo_unpack_dmat(nx, nu, sRSQrq, nu, 0, &H[nx*(nx+nu)], nx+nu);
    blasfeo_unpack_tran_dmat(nx, nu, sRSQrq, nu, 0, &H[nx], nx+nu);

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
    blasfeo_unpack_dvec(nx, srq, nu, &g[0]);
    blasfeo_unpack_dvec(nu, srq, 0, &g[nx]);

    // printf("acados rq (nx = %d, nu = %d)\n", nx, nu);
    // blasfeo_print_tran_dvec(srq->m, srq, 0);
    // printf("qpDUNES g':\n");
    // d_print_mat(1, nx+nu, g, 1);
    // printf("********************************************\n\n");
}



static void form_dynamics(double *ABt, double *b, int nx, int nu, struct blasfeo_dmat *sBAbt, struct blasfeo_dvec *sb)
{
    // copy A
    blasfeo_unpack_dmat(nx, nx, sBAbt, nu, 0, &ABt[0], nx+nu);
    // copy B
    blasfeo_unpack_dmat(nu, nx, sBAbt, 0, 0, &ABt[nx], nx+nu);
    // copy b
    blasfeo_unpack_dvec(nx, sb, 0, b);

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
    for (int ii = 0; ii < nx+nu; ii++)
    {
        zLow[ii] = -infty;
        zUpp[ii] = infty;
    }
    for (int ii = 0; ii < nb; ii++)
    {
        if (idxb[ii] < nu)
        {  // input bound
            zLow[idxb[ii] + nx] = DVECEL_LIBSTR(sd, ii);  // lb[ii]
            zUpp[idxb[ii] + nx] = -DVECEL_LIBSTR(sd, ii + nb + ng);  // ub[ii]
        } else
        {  // state bounds
            zLow[idxb[ii] - nu] = DVECEL_LIBSTR(sd, ii);  // lb[ii]
            zUpp[idxb[ii] - nu] = -DVECEL_LIBSTR(sd, ii + nb + ng);  // ub[ii]
        }
    }
}



static void form_inequalities(double *Ct, double *lc, double *uc, int nx,  int nu, int nb, int ng,
    struct blasfeo_dmat *sDCt, struct blasfeo_dvec *sd)
{
    int ii;

    // copy C
    blasfeo_unpack_dmat(nx, ng, sDCt, nu, 0, &Ct[0], nx+nu);
    // copy D
    blasfeo_unpack_dmat(nu, ng, sDCt, 0, 0, &Ct[nx], nx+nu);
    // copy lc
    blasfeo_unpack_dvec(ng, sd, nb, lc);
    // copy uc
    blasfeo_unpack_dvec(ng, sd, 2*nb+ng, uc);

    for (ii = 0; ii < ng; ii++) uc[ii] = -uc[ii];

}



int ocp_qp_qpdunes_calculate_workspace_size(ocp_qp_dims *dims, void *args_)
{
    int nx = dims->nx[0];
    int nu = dims->nu[0];
    int nDmax = get_maximum_number_of_inequality_constraints(dims);
    int nz = nx + nu;

    int size = sizeof(ocp_qp_qpdunes_workspace);

    size += nz*nz*sizeof(double);  // H
    size += nx*nx*sizeof(double);  // Q
    size += nu*nu*sizeof(double);  // R
    size += nx*nu*sizeof(double);  // S
    size += nz*sizeof(double);  // g
    size += nx*nz*sizeof(double);  // ABt
    size += nz*sizeof(double);  // b
    size += nDmax*nz*sizeof(double);  // Ct
    size += 2*nDmax*sizeof(double);  // lc, uc
    size += 2*nz*sizeof(double);  // zLow, zUpp

    return size;
}



static void cast_workspace(ocp_qp_qpdunes_workspace *work, ocp_qp_qpdunes_memory *mem)
{
    char *c_ptr = (char *)work;
    c_ptr += sizeof(ocp_qp_qpdunes_workspace);

    int nx = mem->nx;
    int nu = mem->nu;
    int nz = mem->nz;
    int nDmax = mem->nDmax;

    assign_double(nz*nz, &work->H, &c_ptr);
    assign_double(nx*nx, &work->Q, &c_ptr);
    assign_double(nu*nu, &work->R, &c_ptr);
    assign_double(nx*nu, &work->S, &c_ptr);
    assign_double(nz, &work->g, &c_ptr);
    assign_double(nx*nz, &work->ABt, &c_ptr);
    assign_double(nz, &work->b, &c_ptr);
    assign_double(nDmax*nz, &work->Ct, &c_ptr);
    assign_double(nDmax, &work->lc, &c_ptr);
    assign_double(nDmax, &work->uc, &c_ptr);
    assign_double(nz, &work->zLow, &c_ptr);
    assign_double(nz, &work->zUpp, &c_ptr);
}



static int update_memory(ocp_qp_in *in, ocp_qp_qpdunes_args *args, ocp_qp_qpdunes_memory *mem,
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
    if (args->warmstart == 0)
        for (int ii = 0; ii < N; ii++)
            for (int jj = 0; jj < nx; jj++)
                mem->qpData.lambda.data[ii*nx+jj] = 0.0;

    mem->qpData.options.maxIter = args->options.maxIter;

    if (mem->firstRun == 1)
    {
        // check if qpDUNES will detect clipping or qpOASES
        stageQps = check_stage_qp_solver(args, in);

        // if user specified clipping but problem requires qpOASES, throw error
        assert(!((args->stageQpSolver == QPDUNES_WITH_CLIPPING)
            && (stageQps == QPDUNES_WITH_QPOASES)) && "Cannot use clipping for this QP");

        // if user specified qpOASES but clipping is detected, trick qpDUNES to detect qpOASES
        // NOTE(dimitris): also needed when partial condensing is used because qpDUNES detects
        // clipping on last stage and crashes..
        if (args->stageQpSolver == QPDUNES_WITH_QPOASES)
        {
            // make Q[N] non diagonal
            DMATEL_LIBSTR(&in->RSQrq[N], 0, 1) += 1e-8;
            DMATEL_LIBSTR(&in->RSQrq[N], 1, 0) += 1e-8;
        }

        // setup of intervals
        for (int kk = 0; kk < N; ++kk) {
            form_RSQ(work->R, work->S, work->Q, nx, nu, &in->RSQrq[kk]);
            form_g(work->g, nx, nu, &in->rq[kk]);
            form_bounds(work->zLow, work->zUpp, nx, nu, nb[kk], ng[kk], in->idxb[kk],
            &in->d[kk], args->options.QPDUNES_INFTY);
            form_dynamics(work->ABt, work->b, nx, nu, &in->BAbt[kk], &in->b[kk]);

            if (args->stageQpSolver == QPDUNES_WITH_QPOASES) {
                if (ng[kk] == 0) {
                    value = qpDUNES_setupRegularInterval(
                        &(mem->qpData), mem->qpData.intervals[kk], 0, work->Q,
                        work->R, work->S, work->g, work->ABt, 0, 0,
                        work->b, work->zLow, work->zUpp, 0, 0, 0, 0, 0, 0, 0);
                } else {
                    form_inequalities(work->Ct, work->lc, work->uc, nx, nu, nb[kk], ng[kk], &in->DCt[kk], &in->d[kk]);
                    value = qpDUNES_setupRegularInterval(
                        &(mem->qpData), mem->qpData.intervals[kk], 0, work->Q,
                        work->R, work->S, work->g, work->ABt, 0, 0,
                        work->b, work->zLow, work->zUpp, 0, 0, 0, 0, work->Ct,
                        work->lc, work->uc);
                }
            } else {  // do not pass S[kk] or Cx[kk]/Cu[kk] at all
                value = qpDUNES_setupRegularInterval(
                    &(mem->qpData), mem->qpData.intervals[kk], 0, work->Q,
                    work->R, 0, work->g, work->ABt, 0, 0, work->b,
                    work->zLow, work->zUpp, 0, 0, 0, 0, 0, 0, 0);
            }
            if (value != QPDUNES_OK) {
                printf("Setup of qpDUNES failed on interval %d\n", kk);
                return (int)value;
            }
        }
        form_bounds(work->zLow, work->zUpp, nx, 0, nb[N], ng[N], in->idxb[N],
            &in->d[N], args->options.QPDUNES_INFTY);
        form_RSQ(work->R, work->S, work->Q, nx, 0, &in->RSQrq[N]);
        form_g(work->g, nx, 0, &in->rq[N]);  // work->g = q
        if (ng[N] == 0) {
            value = qpDUNES_setupFinalInterval(
                &(mem->qpData), mem->qpData.intervals[N], work->Q, work->g,
                work->zLow, work->zUpp, 0, 0, 0);
        } else {
            form_inequalities(work->Ct, work->lc, work->uc, nx, 0, nb[N], ng[N], &in->DCt[N], &in->d[N]);
            value = qpDUNES_setupFinalInterval(
                &(mem->qpData), mem->qpData.intervals[N], work->Q, work->g,
                work->zLow, work->zUpp, work->Ct, work->lc, work->uc);
        }
        if (value != QPDUNES_OK) {
            printf("Setup of qpDUNES failed on last interval\n");
            return (int)value;
        }

        // setup of stage QPs
        value = qpDUNES_setupAllLocalQPs(&(mem->qpData), isLTI = QPDUNES_FALSE);
        if (value != QPDUNES_OK) {
            printf("Setup of qpDUNES failed on initialization of stage QPs\n");
            return (int)value;
        }
    } else {  // if mem->firstRun == 0
        if (args->isLinearMPC == 0) {
            for (int kk = 0; kk < N; kk++) {
                form_H(work->H, nx, nu, &in->RSQrq[kk]);
                form_g(work->g, nx, nu, &in->rq[kk]);
                form_dynamics(work->ABt, work->b, nx, nu, &in->BAbt[kk], &in->b[kk]);

                form_bounds(work->zLow, work->zUpp, nx, nu, nb[kk], ng[kk], in->idxb[kk],
                    &in->d[kk], args->options.QPDUNES_INFTY);
                if (ng[kk] == 0) {
                    value = qpDUNES_updateIntervalData(
                        &(mem->qpData), mem->qpData.intervals[kk], work->H,
                        work->g, work->ABt, work->b, work->zLow, work->zUpp,
                        0, 0, 0, 0);
                } else {
                    form_inequalities(work->Ct, work->lc, work->uc, nx, nu, nb[kk], ng[kk], &in->DCt[kk], &in->d[kk]);
                    value = qpDUNES_updateIntervalData(
                        &(mem->qpData), mem->qpData.intervals[kk], work->H,
                        work->g, work->ABt, work->b, work->zLow, work->zUpp,
                        work->Ct, work->lc, work->uc, 0);
                }
                if (value != QPDUNES_OK) {
                    printf("Update of qpDUNES failed on interval %d\n", kk);
                    return (int)value;
                }
                // qpDUNES_printMatrixData( work->ABt, nx, nx+nu, "AB[%d]", kk);
            }
            form_bounds(work->zLow, work->zUpp, nx, 0, nb[N], ng[N], in->idxb[N],
                &in->d[N], args->options.QPDUNES_INFTY);
            form_RSQ(work->R, work->S, work->Q, nx, 0, &in->RSQrq[N]);
            form_g(work->g, nx, 0, &in->rq[N]);  // work->g = q
            if (ng[N] == 0) {
                value = qpDUNES_updateIntervalData(
                    &(mem->qpData), mem->qpData.intervals[N], work->Q,
                    work->g, 0, 0, work->zLow, work->zUpp, 0, 0, 0, 0);
            } else {
                form_inequalities(work->Ct, work->lc, work->uc, nx, 0, nb[N], ng[N], &in->DCt[N], &in->d[N]);
                value = qpDUNES_updateIntervalData(
                    &(mem->qpData), mem->qpData.intervals[N], work->Q,
                    work->g, 0, 0, work->zLow, work->zUpp, work->Ct, work->lc,
                    work->uc, 0);
            }
            if (value != QPDUNES_OK) {
                printf("Update of qpDUNES failed on last interval\n");
                return (int)value;
            }
        } else {  // linear MPC
            form_bounds(work->zLow, work->zUpp, nx, nu, nb[0], ng[0], in->idxb[0],
                &in->d[0], args->options.QPDUNES_INFTY);
            value = qpDUNES_updateIntervalData(
                &(mem->qpData), mem->qpData.intervals[0], 0, 0, 0, 0,
                work->zLow, work->zUpp, 0, 0, 0, 0);
            if (value != QPDUNES_OK) {
                printf("Update of qpDUNES failed on first interval\n");
                return (int)value;
            }
        }
    }
    mem->firstRun = 0;
    return (int)value;
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
        for (int ii = 0; ii < 2*nv + 2*nc; ii++) dual_sol[ii] = (dual_sol[ii] >= 0.0) ? dual_sol[ii] : 0.0;

        blasfeo_pack_dvec(nx[kk], &mem->qpData.intervals[kk]->z.data[0], &out->ux[kk], nu[kk]);
        blasfeo_pack_dvec(nu[kk], &mem->qpData.intervals[kk]->z.data[nx[kk]], &out->ux[kk], 0);

        for (int ii = 0; ii < 2*nb[kk]+2*ng[kk]; ii++) out->lam[kk].pa[ii] = 0.0;

        // qpdunes multipliers form pairs (lam_low, lam_upp), first for bounds, then for constraints
        for (int ii = 0; ii < nb[kk]; ii++)
        {
            // flip bounds from [x u] to [u x]
            if (ii < nx[kk])
            {
                out->lam[kk].pa[nu[kk]+ii] = dual_sol[2*idxb[kk][ii]];
                out->lam[kk].pa[nb[kk]+ng[kk]+nu[kk]+ii] = dual_sol[2*idxb[kk][ii]+1];
            }
            else
            {
                out->lam[kk].pa[ii-nx[kk]] = dual_sol[2*idxb[kk][ii]];
                out->lam[kk].pa[nb[kk]+ng[kk]+ii-nx[kk]] = dual_sol[2*idxb[kk][ii]+1];
            }

        }

        for (int ii = 0; ii < ng[kk]; ii++)
        {
            out->lam[kk].pa[nb[kk]+ii] = dual_sol[2*nv+2*ii];
            out->lam[kk].pa[2*nb[kk]+ng[kk]+ii] = dual_sol[2*nv+2*ii+1];
        }
    }
    for (int kk = 0; kk < N; kk++)
    {
        blasfeo_pack_dvec(nx[kk+1], &mem->qpData.lambda.data[kk*nx[kk+1]], &out->pi[kk], 0);
    }
}


// TODO(dimitris): free also qp_in, qp_out, etc and write for all solvers?
void ocp_qp_qpdunes_free_memory(void *mem_)
{
    ocp_qp_qpdunes_memory *mem = (ocp_qp_qpdunes_memory *)mem_;
    qpDUNES_cleanup(&(mem->qpData));
}



int ocp_qp_qpdunes(ocp_qp_in *in, ocp_qp_out *out, void *args_, void *mem_, void *work_)
{

    acados_timer tot_timer, qp_timer, interface_timer;
    ocp_qp_info *info = (ocp_qp_info *) out->misc;
    acados_tic(&tot_timer);

    // cast data structures
    ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args *)args_;
    ocp_qp_qpdunes_memory *mem = (ocp_qp_qpdunes_memory *)mem_;
    ocp_qp_qpdunes_workspace *work = (ocp_qp_qpdunes_workspace *)work_;

    acados_tic(&interface_timer);
    cast_workspace(work, mem);
    update_memory(in, args, mem, work);
    info->interface_time = acados_toc(&interface_timer);

    acados_tic(&qp_timer);
    return_t qpdunes_status = qpDUNES_solve(&(mem->qpData));
    info->solve_QP_time = acados_toc(&qp_timer);
    
    acados_tic(&interface_timer);
    fill_in_qp_out(in, out, mem);
    info->interface_time += acados_toc(&interface_timer);

    info->total_time = acados_toc(&tot_timer);
    info->num_iter = mem->qpData.log.numIter;
    
    int acados_status = qpdunes_status;
    if (qpdunes_status == QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND) acados_status = ACADOS_SUCCESS;
    if (qpdunes_status == QPDUNES_ERR_ITERATION_LIMIT_REACHED) acados_status = ACADOS_MAXITER;
    return acados_status;
}
