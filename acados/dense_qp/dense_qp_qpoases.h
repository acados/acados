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


#ifndef ACADOS_DENSE_QP_DENSE_QP_QPOASES_H_
#define ACADOS_DENSE_QP_DENSE_QP_QPOASES_H_

#ifdef __cplusplus
extern "C" {
#endif

// qpoases
#include "qpOASES_e.h"
// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"
// blasfeo
#include "blasfeo_target.h"
#include "blasfeo_common.h"

typedef struct dense_qp_qpoases_args_ {
    double max_cputime;  // maximum cpu time in seconds
    int max_nwsr;        // maximum number of working set recalculations
    int warm_start;      // warm start with dual_sol in memory
} dense_qp_qpoases_args;



typedef struct dense_qp_qpoases_memory_ {
    double *H;
    double *R;
    double *g;
    double *A;
    double *b;
    double *d_lb0;
    double *d_ub0;
    double *d_lb;
    double *d_ub;
    double *C;
    double *d_lg;
    double *d_ug;
    int *idxb;
    double *prim_sol;
    double *dual_sol;
    void *QPB;                 // NOTE(giaf): cast to QProblemB to use
    void *QP;                  // NOTE(giaf): cast to QProblem to use
    double cputime;            // cputime of qpoases
    int nwsr;                  // performed number of working set recalculations
} dense_qp_qpoases_memory;



int dense_qp_qpoases_calculate_args_size(dense_qp_dims *dims);
//
void *dense_qp_qpoases_assign_args(dense_qp_dims *dims, void *raw_memory);
//
void dense_qp_qpoases_initialize_default_args(void *args_);
//
int dense_qp_qpoases_calculate_memory_size(dense_qp_dims *dims, void *args_);
//
void *dense_qp_qpoases_assign_memory(dense_qp_dims *dims, void *args_, void *raw_memory);
//
int dense_qp_qpoases_calculate_workspace_size(dense_qp_dims *dims, void *args_);
//
int dense_qp_qpoases(dense_qp_in *qp_in, dense_qp_out *qp_out, void *args_, void *memory_, void *work_);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_DENSE_QP_DENSE_QP_QPOASES_H_
