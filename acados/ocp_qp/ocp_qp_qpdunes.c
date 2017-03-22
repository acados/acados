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

#include <stdio.h>
#include <stdlib.h>

#include "acados/utils/timing.h"

static int_t max_of_three(int_t a, int_t b, int_t c) {
     int_t ans = a;
     (void)((ans < b) && (ans = b));
     (void)((ans < c) && (ans = c));
     return ans;
}


static void transpose_matrix(real_t *mat, int m, int n, real_t *tmp) {
    int_t ii, jj;

    for (jj = 0; jj < n; jj++) {
        for (ii = 0; ii < m; ii++) {
             tmp[ii*n+jj] = mat[jj*m + ii];
        }
    }
    for (ii = 0; ii < m*n; ii++) mat[ii] = tmp[ii];
}


static void ocp_qp_qpdunes_cast_workspace(ocp_qp_qpdunes_workspace *work,
    ocp_qp_qpdunes_memory *mem) {

    char *ptr = (char *)work;

    ptr += sizeof(ocp_qp_qpdunes_workspace);
    work->At = (real_t*)ptr;
    ptr += (mem->dimA)*sizeof(real_t);
    work->Bt = (real_t*)ptr;
    ptr += (mem->dimB)*sizeof(real_t);
    work->Ct = (real_t*)ptr;
    ptr += (mem->dimC)*sizeof(real_t);
    work->scrap = (real_t*)ptr;
    ptr += (mem->maxDim)*sizeof(real_t);
    work->zLow = (real_t*)ptr;
    ptr += (mem->dimz)*sizeof(real_t);
    work->zUpp = (real_t*)ptr;
    ptr += (mem->dimz)*sizeof(real_t);
    work->g = (real_t*)ptr;
    // ptr += (mem->dimz)*sizeof(real_t);
}


static int_t ocp_qp_qpdunes_update_memory(const ocp_qp_in *in,  const ocp_qp_qpdunes_args *args,
    ocp_qp_qpdunes_memory *mem, ocp_qp_qpdunes_workspace *work) {

    int_t ii, kk, N, nx, nu, nc;
    boolean_t isLTI;  // TODO(dimitris): use isLTI flag for LTI systems
    return_t return_value = 0;

    N = in->N;
    nx = in->nx[0];
    nu = in->nu[0];

    // dummy command
    if (args->options.logLevel == 0) ii = 0;

    if (mem->firstRun == 1) {
        /* setup of intervals */
        for (kk = 0; kk < N; ++kk) {
            for (ii = 0; ii < nx; ii++) work->g[ii] = in->q[kk][ii];
            for (ii = 0; ii < nu; ii++) work->g[ii+nx] = in->r[kk][ii];
            for (ii = 0; ii < nx+nu; ii++) {
                work->zLow[ii] = -args->options.QPDUNES_INFTY;
                work->zUpp[ii] = args->options.QPDUNES_INFTY;
            }
            for (ii = 0; ii < in->nb[kk]; ii++) {
                work->zLow[in->idxb[kk][ii]] = in->lb[kk][ii];
                work->zUpp[in->idxb[kk][ii]] = in->ub[kk][ii];
            }
            for (ii = 0; ii < nx*nx; ii++) work->At[ii] = in->A[kk][ii];
            for (ii = 0; ii < nx*nu; ii++) work->Bt[ii] = in->B[kk][ii];
            transpose_matrix(work->At, nx, nx, work->scrap);
            transpose_matrix(work->Bt, nx, nu, work->scrap);

            if (mem->stageQpSolver == QPDUNES_WITH_QPOASES) {
                nc = in->nc[kk];
                if (nc == 0) {
                    return_value = qpDUNES_setupRegularInterval(&(mem->qpData),
                    mem->qpData.intervals[kk], 0, in->Q[kk], in->R[kk], in->S[kk], work->g, 0,
                    work->At, work->Bt, in->b[kk], work->zLow, work->zUpp, 0, 0, 0, 0, 0, 0, 0);
                } else {
                    for (ii = 0; ii < nc*nx; ii++) work->Ct[ii] = in->Cx[kk][ii];
                    for (ii = 0; ii < nc*nu; ii++) work->Ct[ii+nc*nx] = in->Cu[kk][ii];
                    transpose_matrix(work->Ct, nc, nx+nu, work->scrap);
                    return_value = qpDUNES_setupRegularInterval(&(mem->qpData),
                    mem->qpData.intervals[kk], 0, in->Q[kk], in->R[kk], in->S[kk], work->g, 0,
                    work->At, work->Bt, in->b[kk], work->zLow, work->zUpp, 0, 0, 0, 0,
                    work->Ct, in->lc[kk], in->uc[kk]);
                }
            } else {  // do not pass S[kk] or Cx[kk]/Cu[kk] at all
                return_value = qpDUNES_setupRegularInterval(&(mem->qpData),
                mem->qpData.intervals[kk], 0, in->Q[kk], in->R[kk], 0, work->g, 0,
                work->At, work->Bt, in->b[kk], work->zLow, work->zUpp, 0, 0, 0, 0, 0, 0, 0);
            }
            if (return_value != QPDUNES_OK) {
                printf("Setup of qpDUNES failed on interval %d\n", kk);
                return (int_t)return_value;
            }
        }
        for (ii = 0; ii < nx; ii++) {
            work->zLow[ii] = -args->options.QPDUNES_INFTY;
            work->zUpp[ii] = args->options.QPDUNES_INFTY;
        }
        for (ii = 0; ii < in->nb[N]; ii++) {
            work->zLow[in->idxb[N][ii]] = in->lb[N][ii];
            work->zUpp[in->idxb[N][ii]] = in->ub[N][ii];
        }
        nc = in->nc[N];
        if (nc == 0) {
            return_value = qpDUNES_setupFinalInterval(&(mem->qpData), mem->qpData.intervals[N],
            in->Q[N], in->q[N], work->zLow, work->zUpp, 0, 0, 0);
        } else {
            for (ii = 0; ii < nc*nx; ii++) work->Ct[ii] = in->Cx[N][ii];
            transpose_matrix(work->Ct, nc, nx, work->scrap);
            return_value = qpDUNES_setupFinalInterval(&(mem->qpData), mem->qpData.intervals[N],
            in->Q[N], in->q[N], work->zLow, work->zUpp, work->Ct, in->lc[N], in->uc[N]);
        }

        if (return_value != QPDUNES_OK) {
            printf("Setup of qpDUNES failed on last interval\n");
            return (int_t)return_value;
        }

        /* setup of stage QPs */
        return_value = qpDUNES_setupAllLocalQPs(&(mem->qpData), isLTI = QPDUNES_FALSE);
        if (return_value != QPDUNES_OK) {
            printf("Setup of qpDUNES failed on initialization of stage QPs\n");
            return (int_t)return_value;
        }
    } else {  // if mem->firstRun == 0
        // TODO(dimitris): update data if NOT first run!
    }

    mem->firstRun = 0;
    return (int_t)return_value;
}


static void fill_in_qp_out(const ocp_qp_in *in, ocp_qp_out *out, ocp_qp_qpdunes_memory *mem) {
    int ii, kk;

    for (kk = 0; kk < in->N+1; kk++) {
        for (ii = 0; ii < in->nx[kk]; ii++) {
            out->x[kk][ii] = mem->qpData.intervals[kk]->z.data[ii];
        }
        for (ii = 0; ii < in->nu[kk]; ii++) {
            out->u[kk][ii] = mem->qpData.intervals[kk]->z.data[in->nx[kk]+ii];
        }
    }
    // TODO(dimitris): fill-in multipliers
}


static qpdunes_stage_qp_solver_t define_stage_qp_solver(const ocp_qp_in *in) {
    int_t ii, jj, kk;
    int_t nD = 0;

    // check for polyhedral constraints
    for (kk = 0; kk < in->N+1; kk++) nD += in->nc[kk];
    if (nD != 0) return QPDUNES_WITH_QPOASES;

    // check for non-zero cross terms
    for (kk = 0; kk < in->N; kk++) {
        for (ii = 0; ii < in->nx[kk]*in->nu[kk]; ii++) {
            if (in->S[kk][ii] != 0) {
                return QPDUNES_WITH_QPOASES;
            }
        }
    }

    // check for non-diagonal Q
    for (kk = 0; kk < in->N+1; kk++) {
        for (jj = 0; jj < in->nx[kk]; jj++) {
            for (ii = 0; ii < in->nx[kk]; ii++) {
                if ((ii != jj) && (in->Q[kk][jj*in->nx[kk]+ii] != 0)) {
                    return QPDUNES_WITH_QPOASES;
                }
            }
        }
    }

    // check for non-diagonal R
    for (kk = 0; kk < in->N; kk++) {
        for (jj = 0; jj < in->nu[kk]; jj++) {
            for (ii = 0; ii < in->nu[kk]; ii++) {
                if ((ii != jj) && (in->R[kk][jj*in->nu[kk]+ii] != 0)) {
                    return QPDUNES_WITH_QPOASES;
                }
            }
        }
    }

    return QPDUNES_WITH_CLIPPING;
}


static int_t get_maximum_number_of_inequality_constraints(const ocp_qp_in *in) {
    int_t kk, nDmax;

    nDmax = in->nc[0];
    for (kk = 1; kk < in->N+1; kk ++) {
        if (in->nc[kk] > nDmax) nDmax = in->nc[kk];
    }
    return nDmax;
}


int_t ocp_qp_qpdunes_create_arguments(void *args_, int_t opts_) {
    ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args*) args_;
    qpdunes_options_t opts = (qpdunes_options_t) opts_;

    if (opts == QPDUNES_DEFAULT_ARGUMENTS) {
        args->options = qpDUNES_setupDefaultOptions();
        args->options.printLevel = 0;
    } else {
        printf("\nUknown option (%d) for qpDUNES!\n", opts_);
        return -1;
    }
    return 0;
}


int_t ocp_qp_qpdunes_calculate_workspace_size(const ocp_qp_in *in, void *args_) {
    ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args*) args_;

    int_t size, dimA, dimB, dimC, nDmax, dimz, maxDim;

    // dummy command
    if (args->options.logLevel == 0) dimA = 0;

    nDmax = get_maximum_number_of_inequality_constraints(in);
    dimA = in->nx[0]*in->nx[0];
    dimB = in->nx[0]*in->nu[0];
    dimz = in->nx[0]+in->nu[0];
    dimC = nDmax*dimz;

    // calculate memory size for scrap memory (used to transpose matrices)
    maxDim = max_of_three(dimA, dimB, dimC);

    size = sizeof(ocp_qp_qpdunes_workspace);
    size += (dimA + dimB + dimC + maxDim + 3*dimz)*sizeof(real_t);
    return size;
}


int_t ocp_qp_qpdunes_create_memory(const ocp_qp_in *in, void *args_, void *mem_) {
    ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args*) args_;
    ocp_qp_qpdunes_memory *mem = (ocp_qp_qpdunes_memory *) mem_;

    int_t N, nx, nu, kk;
    uint_t *nD_ptr = 0;
    return_t return_value;

    N = in->N;
    nx = in->nx[0];
    nu = in->nu[0];

    mem->firstRun = 1;
    mem->dimA = nx*nx;
    mem->dimB = nx*nu;
    mem->dimz = nx+nu;
    mem->nDmax = get_maximum_number_of_inequality_constraints(in);
    mem->dimC = mem->nDmax*mem->dimz;
    mem->maxDim = max_of_three(mem->dimA, mem->dimB, mem->dimC);

    /* Check for constant dimensions */
    for (kk = 1; kk < N; kk++) {
        if ((nx != in->nx[kk]) || (nu != in->nu[kk])) {
            printf("\nqpDUNES does not support varying dimensions!");
            return -1;
        }
    }
    if ((nx != in->nx[N]) || (in->nu[N] != 0)) {
        printf("\nqpDUNES does not support varying dimensions!");
        return -1;
    }

    mem->stageQpSolver = define_stage_qp_solver(in);
    // if (mem->stageQpSolver == QPDUNES_WITH_QPOASES) printf("\n\n >>>>>> QPDUNES + QPOASES!\n");
    // if (mem->stageQpSolver == QPDUNES_WITH_CLIPPING) printf("\n\n >>>>>> QPDUNES + CLIPPING!\n");

    if (mem->stageQpSolver == QPDUNES_WITH_QPOASES) {
        // NOTE: lsType 5 seems to work but yields wrong results with ineq. constraints
        if (args->options.lsType != 7) {
            args->options.lsType = 7;
            // TODO(dimitris): write proper acados warnings and errors
            if (args->options.printLevel > 0) {
                printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                printf("WARNING: Changed line-search algorithm for qpDUNES (incompatible with QP)");
                printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            }
        }
        if (mem->nDmax > 0) nD_ptr = (uint_t*)in->nc;  // otherwise leave pointer equal to zero
    }

    /* memory allocation */
    return_value = qpDUNES_setup(&(mem->qpData), N, nx, nu, nD_ptr, &(args->options));
    if (return_value != QPDUNES_OK) {
        printf("Setup of the QP solver failed\n");
        return (int_t)return_value;
    }
    return 0;
}


void ocp_qp_qpdunes_free_memory(void *mem_) {
    ocp_qp_qpdunes_memory *mem = (ocp_qp_qpdunes_memory *) mem_;
    qpDUNES_cleanup(&(mem->qpData));
}


int_t ocp_qp_qpdunes(ocp_qp_in *in, ocp_qp_out *out, void *args_, void *mem_, void *work_) {
    ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args*) args_;
    ocp_qp_qpdunes_memory *mem = (ocp_qp_qpdunes_memory *) mem_;
    ocp_qp_qpdunes_workspace *work = (ocp_qp_qpdunes_workspace *) work_;

    return_t return_value;

    ocp_qp_qpdunes_cast_workspace(work, mem);
    ocp_qp_qpdunes_update_memory(in, args, mem, work);

    return_value = qpDUNES_solve(&(mem->qpData));
    if (return_value != QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND) {
        printf("qpDUNES failed to solve the QP. Error code: %d\n", return_value);
        return (int_t)return_value;
    }
    fill_in_qp_out(in, out, mem);

    return 0;
}
