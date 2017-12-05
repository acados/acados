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


#ifndef ACADOS_DENSE_QP_DENSE_QP_QORE_H_
#define ACADOS_DENSE_QP_DENSE_QP_QORE_H_

#ifdef __cplusplus
extern "C" {
#endif

// qore
#include "qpsolver_dense.h"
// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"

typedef struct dense_qp_qore_args_ {
    int nsmax;          // maximum size of Schur complement
    int prtfreq;        // print frequency,
                        // prtfreq  < 0: disable printing;
                        // prtfreq == 0: print on each call and include working set changes;
                        // prtfreq  > 0: print on every prtfreq seconds, but do not include working set changes;
    int warm_start;     // warm start with updated matrices H and C
    int warm_strategy;  // 0: ramp-up from zero homotopy; 1: setup homotopy from the previous solution
    int hot_start;      // hot start with unchanged matrices H and C
} dense_qp_qore_args;



typedef struct dense_qp_qore_memory_ {
    double *H;
    double *g;
    double *A;
    double *b;
    double *C;
    double *Ct;
    double *d_lb0;
    double *d_ub0;
    double *d_lb;
    double *d_ub;
    double *d_lg;
    double *d_ug;
    double *lb;
    double *ub;
    int *idxb;
    double *prim_sol;
    double *dual_sol;
    QoreProblemDense *QP;
    int num_iter;
} dense_qp_qore_memory;



int dense_qp_qore_calculate_args_size(dense_qp_dims *dims);
//
void *dense_qp_qore_assign_args(dense_qp_dims *dims, void *raw_memory);
//
void dense_qp_qore_initialize_default_args(void *args_);
//
int dense_qp_qore_calculate_memory_size(dense_qp_dims *dims, void *args_);
//
void *dense_qp_qore_assign_memory(dense_qp_dims *dims, void *args_, void *raw_memory);
//
int dense_qp_qore_calculate_workspace_size(dense_qp_dims *dims, void *args_);
//
int dense_qp_qore(dense_qp_in *qp_in, dense_qp_out *qp_out, void *args_, void *memory_, void *work_);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_DENSE_QP_DENSE_QP_QORE_H_
