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

#include <stdio.h>
#include <stdlib.h>
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_i_aux.h"

#include "qpDUNES/include/qpDUNES.h"
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#include "acados/utils/timing.h"

#define qpDUNES_INFTY 1e10

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
    int_t ii, kk, N, nx, nu;
    boolean_t isLTI;
    return_t return_value;

    N = in->N;
    nx = in->nx[0];
    nu = in->nu[0];

    // dummy command
    if (args->options.logLevel == 0) ii = 0;

    if (mem->firstRun == 1) {
        /* setup of intervals */
        for (kk = 0; kk < N; ++kk) {
            // TODO(dimitris): move to a function
            for (ii = 0; ii < nx; ii++) work->g[ii] = in->q[kk][ii];
            for (ii = 0; ii < nu; ii++) work->g[ii+nx] = in->r[kk][ii];
            for (ii = 0; ii < nx+nu; ii++) {
                work->zLow[ii] = -qpDUNES_INFTY;
                work->zUpp[ii] = qpDUNES_INFTY;
            }
            for (ii = 0; ii < in->nb[kk]; ii++) {
                // TODO(dimitris): what's our infty?
                work->zLow[in->idxb[kk][ii]] = in->lb[kk][ii];
                work->zUpp[in->idxb[kk][ii]] = in->ub[kk][ii];
            }
            for (ii = 0; ii < nx*nx; ii++) work->At[ii] = in->A[kk][ii];
            for (ii = 0; ii < nx*nu; ii++) work->Bt[ii] = in->B[kk][ii];
            transpose_matrix(work->At, nx, nx, work->scrap);
            transpose_matrix(work->Bt, nx, nu, work->scrap);
            return_value = qpDUNES_setupRegularInterval(&(mem->qpData), mem->qpData.intervals[kk],
            0, in->Q[kk], in->R[kk], 0, work->g, 0, work->At, work->Bt, in->b[kk], work->zLow,
            work->zUpp, 0, 0, 0, 0, 0, 0, 0);

            if (return_value != QPDUNES_OK) {
                printf("Setup of qpDUNES failed on interval %d\n", kk);
                return (int_t)return_value;
            }
        }
        // TODO(dimitris): can I merge this above?
        for (ii = 0; ii < nx; ii++) {
            work->zLow[ii] = -qpDUNES_INFTY;
            work->zUpp[ii] = qpDUNES_INFTY;
        }
        for (ii = 0; ii < in->nb[N]; ii++) {
            work->zLow[in->idxb[N][ii]] = in->lb[N][ii];
            work->zUpp[in->idxb[N][ii]] = in->ub[N][ii];
        }
        return_value = qpDUNES_setupFinalInterval(&(mem->qpData), mem->qpData.intervals[N],
        in->Q[N], in->q[N], work->zLow, work->zUpp, 0, 0, 0);

        if (return_value != QPDUNES_OK) {
            printf("Setup of qpDUNES failed on last interval\n");
            return (int_t)return_value;
        }

        /* setup of stage QPs */
        // TODO(dimitris): check isLTI flag
        return_value = qpDUNES_setupAllLocalQPs(&(mem->qpData), isLTI = QPDUNES_FALSE);
        if (return_value != QPDUNES_OK) {
            printf("Setup of qpDUNES failed on initialization of stage QPs\n");
            return (int_t)return_value;
        }
    }
    mem->firstRun = 0;
    return (int_t)return_value;
}


static void fill_in_qp_out(ocp_qp_in *in, ocp_qp_out *out, ocp_qp_qpdunes_memory *mem) {
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


int_t ocp_qp_qpdunes_create_arguments(void *args_, int_t opts_) {
    ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args*) args_;
    qpdunes_options_t opts = (qpdunes_options_t) opts_;

    if (opts == QPDUNES_DEFAULT_ARGUMENTS) {
        args->options = qpDUNES_setupDefaultOptions();
    } else {
        printf("\nUknown option (%d) for qpDUNES!\n", opts_);
        return -1;
    }
    return 0;
}


int_t ocp_qp_qpdunes_calculate_workspace_size(const ocp_qp_in *in, void *args_) {
    ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args*) args_;

    int_t size, dimA, dimB, dimz, maxDim;

    // dummy command
    if (args->options.logLevel == 0) dimA = 0;

    dimA = in->nx[0]*in->nx[0];
    dimB = in->nx[0]*in->nu[0];
    dimz = in->nx[0] + in->nu[0];

    maxDim = dimA;
    if (dimB > maxDim) maxDim = dimB;

    size = sizeof(ocp_qp_qpdunes_workspace);
    size += (dimA + dimB + maxDim + 3*dimz)*sizeof(real_t);
    return size;
}


int_t ocp_qp_qpdunes_create_memory(const ocp_qp_in *in, void *args_, void *mem_) {
    ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args*) args_;
    ocp_qp_qpdunes_memory *mem = (ocp_qp_qpdunes_memory *) mem_;
    int_t N, nx, nu, kk;
    uint_t *nD_ptr = 0;
    int_t nD = 0;  // number of ineq. constraints (NOTE: has to be zero for clipping)
    return_t return_value;

    N = in->N;
    nx = in->nx[0];
    nu = in->nu[0];

    mem->firstRun = 1;
    mem->dimA = nx*nx;
    mem->dimB = nx*nu;
    mem->dimz = nx+nu;
    mem->maxDim = mem->dimA;
    if (mem->dimB > mem->maxDim) mem->maxDim = mem->dimB;

    /* Check for constant dimensions and conditions for qpDUNES+clipping */
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
    for (kk = 0; kk < N+1; kk++) nD += in->nc[kk];
    if (nD != 0) {
        printf("qpDUNES with clipping does not support inequality constraints!\n");
        return -1;
    }

    /* memory allocation */
    return_value = qpDUNES_setup(&(mem->qpData), N, nx, nu, nD_ptr, &(args->options));
    if (return_value != QPDUNES_OK) {
        printf("Setup of the QP solver failed\n");
        return (int_t)return_value;
    }
    return (int_t)return_value;
}


void ocp_qp_qpdunes_free_memory(void *mem_) {
    ocp_qp_qpdunes_memory *mem = (ocp_qp_qpdunes_memory *) mem_;
    qpDUNES_cleanup(&(mem->qpData));
}


// TODO(dimitris): Move casts after var declarations here and in OOQP
int_t ocp_qp_qpdunes(ocp_qp_in *in, ocp_qp_out *out, void *args_, void *mem_, void *work_) {
    return_t return_value;

    ocp_qp_qpdunes_args *args = (ocp_qp_qpdunes_args*) args_;
    ocp_qp_qpdunes_memory *mem = (ocp_qp_qpdunes_memory *) mem_;
    ocp_qp_qpdunes_workspace *work = (ocp_qp_qpdunes_workspace *) work_;

    ocp_qp_qpdunes_cast_workspace(work, mem);
    ocp_qp_qpdunes_update_memory(in, args, mem, work);

    // dummy commands
    if (mem->firstRun || args->options.logLevel == 1) work->tmp = 31;
    if (in->nx[0] == 1) out->x[0][0] = 1;

    return_value = qpDUNES_solve(&(mem->qpData));
    if (return_value != QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND) {
        printf("qpDUNES failed to solve the QP. Error code: %d\n", return_value);
        return (int_t)return_value;
    }
    fill_in_qp_out(in, out, mem);

    return return_value;
}
