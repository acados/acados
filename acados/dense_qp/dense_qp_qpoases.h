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

// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"

typedef struct dense_qp_qpoases_opts_
{
    double max_cputime;  // maximum cpu time in seconds
    int max_nwsr;        // maximum number of working set recalculations
    int warm_start;      // warm start with dual_sol in memory
    int use_precomputed_cholesky;
    int hotstart;  // this option requires constant data matrices! (eg linear MPC, inexact schemes
                   // with frozen sensitivities)
    int set_acado_opts;  // use same options as in acado code generation
    int compute_t;       // compute t in qp_out (to have correct residuals in NLP)
} dense_qp_qpoases_opts;

typedef struct dense_qp_qpoases_memory_
{
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
    void *QPB;       // NOTE(giaf): cast to QProblemB to use
    void *QP;        // NOTE(giaf): cast to QProblem to use
    double cputime;  // cputime of qpoases
    int nwsr;        // performed number of working set recalculations
    int first_it;    // to be used with hotstart
} dense_qp_qpoases_memory;

int dense_qp_qpoases_opts_calculate_size(void *config, dense_qp_dims *dims);
//
void *dense_qp_qpoases_opts_assign(void *config, dense_qp_dims *dims, void *raw_memory);
//
void dense_qp_qpoases_opts_initialize_default(void *config, dense_qp_dims *dims, void *opts_);
//
void dense_qp_qpoases_opts_update(void *config, dense_qp_dims *dims, void *opts_);
//
int dense_qp_qpoases__memorycalculate_size(void *config, dense_qp_dims *dims, void *opts_);
//
void *dense_qp_qpoases__memoryassign(void *config, dense_qp_dims *dims, void *opts_,
                                     void *raw_memory);
//
int dense_qp_qpoases__workspacecalculate_size(void *config, dense_qp_dims *dims, void *opts_);
//
int dense_qp_qpoases(void *config, dense_qp_in *qp_in, dense_qp_out *qp_out, void *opts_,
                     void *memory_, void *work_);
//
void dense_qp_qpoases_config_initialize_default(void *config_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_DENSE_QP_DENSE_QP_QPOASES_H_
