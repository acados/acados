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



// static int max_of_two(int a, int b)
// {
//     int ans = a;
//     (void)((ans < b) && (ans = b));
//     return ans;
// }



static int get_maximum_number_of_inequality_constraints(ocp_qp_dims *dims)
{
    int nDmax = dims->ng[0];

    for (int kk = 1; kk < dims->N + 1; kk++)
    {
        if (dims->ng[kk] > nDmax) nDmax = dims->ng[kk];
    }
    return nDmax;
}



static bool check_stage_qp_solver(ocp_qp_qpdunes_args *args, ocp_qp_in *qp_in)
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
    for (int kk = 0; kk < N+1; kk++)
    {
        for (int ii = 0; ii < nx; ii++)
        {
            for (int jj = 0; jj < nx; jj++)
            {
                if ((ii != jj) && (DMATEL_LIBSTR(&qp_in->RSQrq[kk], ii+nu, jj+nu) != 0))
                {
                    stageQpSolver = QPDUNES_WITH_QPOASES;
                }
            }
        }
    }

    // check for non-diagonal R
    for (int kk = 0; kk < N+1; kk++)
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
    if (stageQpSolver == QPDUNES_WITH_QPOASES)
    {
        // NOTE(dimitris): the opposite is possible (i.e., use qpOASES for a clipping formulation)
        if (args->stageQpSolver == QPDUNES_WITH_QPOASES)
        {
            return true;
        } else
        {
            return false;
        }
    }
    return true;
}



int ocp_qp_qpdunes_calculate_args_size(ocp_qp_dims *dims)
{
    int size = 0;
    size += sizeof(ocp_qp_qpdunes_args);
    return size;
}



void *ocp_qp_qpdunes_assign_args(ocp_qp_dims *dims, void *raw_memory)
{
    ocp_qp_qpdunes_args *args;

    char *c_ptr = (char *) raw_memory;

    args = (ocp_qp_qpdunes_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_qpdunes_args);

    assert((char*)raw_memory + ocp_qp_qpdunes_calculate_args_size(dims) == c_ptr);

    return (void *)args;
}



void ocp_qp_qpdunes_initialize_default_args(void *args_)
{
    ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args *)args_;

    qpdunes_options_t opts = QPDUNES_NONLINEAR_MPC;

    args->stageQpSolver = QPDUNES_WITH_CLIPPING;

    if (opts == QPDUNES_DEFAULT_ARGUMENTS) {
        args->options = qpDUNES_setupDefaultOptions();
        args->isLinearMPC = 0;
    } else if (opts == QPDUNES_NONLINEAR_MPC) {
        args->options = qpDUNES_setupDefaultOptions();
        args->isLinearMPC = 0;
        args->options.printLevel = 0;
    } else if (opts == QPDUNES_LINEAR_MPC) {
        args->options = qpDUNES_setupDefaultOptions();
        args->isLinearMPC = 1;
        args->options.printLevel = 0;
    } else {
        printf("\nUnknown option (%d) for qpDUNES!\n", opts);
        args->options = qpDUNES_setupDefaultOptions();
        args->isLinearMPC = 0;
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
    mem->dimA = nx * nx;
    mem->dimB = nx * nu;
    mem->dimz = nx + nu;
    mem->dimC = mem->nDmax * mem->dimz;
    // mem->maxDim = max_of_two(mem->dimA + mem->dimB, mem->dimC);

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


// static void transpose_matrix(double *mat, int m, int n, double *tmp)
// {
//     for (int jj = 0; jj < n; jj++) {
//         for (int ii = 0; ii < m; ii++) {
//             tmp[ii * n + jj] = mat[jj * m + ii];
//         }
//     }
//     for (int ii = 0; ii < m * n; ii++) mat[ii] = tmp[ii];
// }


static void form_H(double *H, int nx, int nu, struct d_strmat *sRSQrq)
{
    // copy Q
    d_cvt_strmat2mat(nx, nx, sRSQrq, nu, nu, &H[0], nx+nu);
    // copy R
    d_cvt_strmat2mat(nu, nu, sRSQrq, 0, 0, &H[nx*(nx+nu) + nx], nx+nu);
    // copy S
    // TODO(dimitris): NOT SURE ABOUT DIMENSIONS OF S AND WHERE TO TRANSPOSE!
    d_cvt_strmat2mat(nu, nx, sRSQrq, nu, 0, &H[nx], nx+nu);
    d_cvt_tran_strmat2mat(nu, nx, sRSQrq, nu, 0, &H[nx*nx + nx*nu], nx+nu);

    // printf("Hessian:\n");
    // d_print_mat(nx+nu, nx+nu, H, nx+nu);
    // exit(1);
}



static void form_RSQ(double *R, double *S, double *Q, int nx, int nu, struct d_strmat *sRSQrq)
{
    // copy Q
    d_cvt_strmat2mat(nx, nx, sRSQrq, nu, nu, Q, nx);
    // copy R
    d_cvt_strmat2mat(nu, nu, sRSQrq, 0, 0, R, nu);
    // copy S
    // TODO(dimitris): NOT SURE ABOUT DIMENSIONS OF S AND WHERE TO TRANSPOSE!
    d_cvt_strmat2mat(nu, nx, sRSQrq, nu, 0, S, nu);
}




static void form_g(double *g, int nx, int nu, struct d_strvec *srq)
{
    d_cvt_strvec2vec(nx, srq, nu, &g[0]);
    d_cvt_strvec2vec(nu, srq, 0, &g[nx]);
}



static void form_dynamics(double *ABt, double *b, int nx, int nu, struct d_strmat *sBAbt, struct d_strvec *sb)
{
    d_cvt_strmat2mat(nx, nx, sBAbt, 0, nu, &ABt[0], nx);
    d_cvt_strmat2mat(nx, nu, sBAbt, 0, 0, &ABt[nx*nx], nx);
    d_cvt_strvec2vec(nx, sb, 0, b);
}



static void form_bounds(double *zLow, double *zUpp, int nx, int nu, int nb, int ng, int *idxb,
    struct d_strvec *sd, double infty)
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
            zUpp[idxb[ii] + nx] = DVECEL_LIBSTR(sd, ii + nb + ng);  // ub[ii]
        } else
        {  // state bounds
            zLow[idxb[ii] - nu] = DVECEL_LIBSTR(sd, ii);  // lbk[ii]
            zUpp[idxb[ii] - nu] = DVECEL_LIBSTR(sd, ii + nb + ng);  // ub[ii]
        }
    }
}



static void form_inequalities(double *Ct, double *lc, double *uc, int nx,  int nu, int nb, int ng,
    struct d_strmat *sDCt, struct d_strvec *sd)
{
    // copy C
    d_cvt_strmat2mat(nx, ng, sDCt, nu, 0, &Ct[0], nx+nu);
    // copy D
    d_cvt_strmat2mat(nu, ng, sDCt, 0, 0, &Ct[nx*ng], nx+nu);
    // copy lc
    d_cvt_strvec2vec(ng, sd, nb, lc);
    // copy uc
    d_cvt_strvec2vec(ng, sd, 2*nb+ng, uc);

}



int ocp_qp_qpdunes_calculate_workspace_size(ocp_qp_dims *dims, void *args_)
{
    ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args *)args_;

    int N = dims->N;
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

    assert(check_stage_qp_solver(args, in) == true);

    int N = in->dim->N;
    int nx = in->dim->nx[0];
    int nu = in->dim->nu[0];
    int *nb = in->dim->nb;
    int *ng = in->dim->ng;

    if (mem->firstRun == 1) {
        /* setup of intervals */
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

        /* setup of stage QPs */
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
                // qpDUNES_printMatrixData( work->ABt, nx, nx+nu, "AB[%d]", kk
                // );
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



static void fill_in_qp_out(ocp_qp_dims *dims, ocp_qp_out *out, ocp_qp_qpdunes_memory *mem)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;

    for (int kk = 0; kk < N + 1; kk++) {
        d_cvt_vec2strvec(nx[kk], &mem->qpData.intervals[kk]->z.data[0], &out->ux[kk], nu[kk]);
        d_cvt_vec2strvec(nu[kk], &mem->qpData.intervals[kk]->z.data[nx[kk]], &out->ux[kk], 0);
        // for (int ii = 0; ii < nx[kk]; ii++) {
        //     out->x[kk][ii] = mem->qpData.intervals[kk]->z.data[ii];
        // }
        // for (int ii = 0; ii < nu[kk]; ii++) {
        //     out->u[kk][ii] = mem->qpData.intervals[kk]->z.data[nx[kk] + ii];
        // }
    }
    // int nn = 0;
    for (int kk = 0; kk < N; kk++) {
        d_cvt_vec2strvec(nx[kk+1], mem->qpData.lambda.data, &out->pi[kk], 0);
        // for (int ii = 0; ii < nx[kk + 1]; ii++) {
        //     out->pi[kk][ii] = mem->qpData.lambda.data[nn++];
        // }
    }
    // TODO(dimitris): fill-in multipliers for inequalities
}


// ocp_qp_qpdunes_args *ocp_qp_qpdunes_create_arguments(qpdunes_options_t opts) {
//     ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args *) malloc(sizeof(ocp_qp_qpdunes_args));

//     if (opts == QPDUNES_DEFAULT_ARGUMENTS) {
//         args->options = qpDUNES_setupDefaultOptions();
//         args->isLinearMPC = 0;
//     } else if (opts == QPDUNES_NONLINEAR_MPC) {
//         args->options = qpDUNES_setupDefaultOptions();
//         args->isLinearMPC = 0;
//         args->options.printLevel = 0;
//     } else if (opts == QPDUNES_LINEAR_MPC) {
//         args->options = qpDUNES_setupDefaultOptions();
//         args->isLinearMPC = 1;
//         args->options.printLevel = 0;
//     } else {
//         printf("\nUnknown option (%d) for qpDUNES!\n", opts);
//         args->options = qpDUNES_setupDefaultOptions();
//         args->isLinearMPC = 0;
//     }
//     return args;
// }


// ocp_qp_qpdunes_memory *ocp_qp_qpdunes_create_memory(const ocp_qp_in *in, void *args_) {

//     ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args *) args_;
//     ocp_qp_qpdunes_memory *mem = (ocp_qp_qpdunes_memory *) malloc(sizeof(ocp_qp_qpdunes_memory));
//     int N, nx, nu, kk;
//     uint *nD_ptr = 0;
//     return_t return_value;

//     N = in->N;
//     nx = in->nx[0];
//     nu = in->nu[0];

//     mem->firstRun = 1;
//     mem->dimA = nx * nx;
//     mem->dimB = nx * nu;
//     mem->dimz = nx + nu;
//     mem->nDmax = get_maximum_number_of_inequality_constraints(in);
//     mem->dimC = mem->nDmax * mem->dimz;
//     mem->maxDim = max_of_two(mem->dimA + mem->dimB, mem->dimC);

//     /* Check for constant dimensions */
//     for (kk = 1; kk < N; kk++) {
//         if ((nx != in->nx[kk]) || (nu != in->nu[kk])) {
//             printf("\nqpDUNES does not support varying dimensions!\n");
//             free(mem);
//             return NULL;
//         }
//     }
//     if ((nx != in->nx[N]) || (in->nu[N] != 0)) {
//         printf("\nqpDUNES does not support varying dimensions!\n");
//         free(mem);
//         return NULL;
//     }

//     mem->stageQpSolver = define_stage_qp_solver(in);
//     // if (mem->stageQpSolver == QPDUNES_WITH_QPOASES) printf("\n\n >>>>>>
//     // QPDUNES + QPOASES!\n"); if (mem->stageQpSolver == QPDUNES_WITH_CLIPPING)
//     // printf("\n\n >>>>>> QPDUNES + CLIPPING!\n");

//     if (mem->stageQpSolver == QPDUNES_WITH_QPOASES) {
//         // NOTE: lsType 5 seems to work but yields wrong results with ineq.
//         // constraints
//         if (args->options.lsType != 7) {
//             args->options.lsType = 7;
//             // TODO(dimitris): write proper acados warnings and errors
//             if (args->options.printLevel > 0) {
//                 printf(
//                     "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
//                     "!!!!!!!!!!!!!\n");
//                 printf(
//                     "WARNING: Changed line-search algorithm for qpDUNES "
//                     "(incompatible with QP)");
//                 printf(
//                     "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
//                     "!!!!!!!!!!!!!\n");
//             }
//         }
//         if (mem->nDmax > 0)
//             nD_ptr = (uint *)in->nc;  // otherwise leave pointer equal to zero
//     }

//     /* memory allocation */
//     return_value = qpDUNES_setup(&(mem->qpData), N, nx, nu, nD_ptr, &(args->options));
//     if (return_value != QPDUNES_OK) {
//         printf("Setup of the QP solver failed\n");
//         free(mem);
//         return NULL;
//     }
//     return mem;
// }



void ocp_qp_qpdunes_free_memory(void *mem_)
{
    ocp_qp_qpdunes_memory *mem = (ocp_qp_qpdunes_memory *)mem_;
    qpDUNES_cleanup(&(mem->qpData));
}



int ocp_qp_qpdunes(ocp_qp_in *in, ocp_qp_out *out, void *args_, void *mem_, void *work_)
{
    ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args *)args_;
    ocp_qp_qpdunes_memory *mem = (ocp_qp_qpdunes_memory *)mem_;
    ocp_qp_qpdunes_workspace *work = (ocp_qp_qpdunes_workspace *)work_;

    return_t return_value;
    // printf("$$ FIRST RUN FLAG %d\n", mem->firstRun);

    cast_workspace(work, mem);

    update_memory(in, args, mem, work);

    return_value = qpDUNES_solve(&(mem->qpData));

    if (return_value != QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND) {
        printf("qpDUNES failed to solve the QP. Error code: %d\n", return_value);
        return (int)return_value;
    }

    fill_in_qp_out(in->dim, out, mem);

    // TODO(dimitris): use acados return value
    return 0;
}