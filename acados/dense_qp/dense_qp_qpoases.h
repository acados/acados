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

/* Ignore compiler warnings from qpOASES */
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wtautological-pointer-compare"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wunused-function"
#include "qpOASES_e/QProblemB.h"
#include "qpOASES_e/QProblem.h"
#pragma clang diagnostic pop
#elif defined(__GNUC__)
    #if __GNUC__ >= 6
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
        #pragma GCC diagnostic ignored "-Wunused-parameter"
        #pragma GCC diagnostic ignored "-Wunused-function"
        #include "qpOASES_e/QProblemB.h"
        #include "qpOASES_e/QProblem.h"
        #pragma GCC diagnostic pop
    #else
        #pragma GCC diagnostic ignored "-Wunused-parameter"
        #pragma GCC diagnostic ignored "-Wunused-function"
        #include "qpOASES_e/QProblemB.h"
        #include "qpOASES_e/QProblem.h"
    #endif
#else
    #include "qpOASES_e/QProblemB.h"
    #include "qpOASES_e/QProblem.h"
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
	int use_precomputed_cholesky;
	int hotstart; 		 // this option requires constant data matrices! (eg linear MPC, inexact schemes with frozen sensitivities) 
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
	int first_it;              // to be used with hotstart
} dense_qp_qpoases_memory;



int dense_qp_qpoases_calculate_args_size(dense_qp_dims *dims, void *submodules_);
//
void *dense_qp_qpoases_assign_args(dense_qp_dims *dims, void **submodules_, void *raw_memory);
//
void *dense_qp_qpoases_copy_args(dense_qp_dims *dims, void *raw_memory, void *source_);
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
