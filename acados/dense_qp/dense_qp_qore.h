/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef ACADOS_DENSE_QP_DENSE_QP_QORE_H_
#define ACADOS_DENSE_QP_DENSE_QP_QORE_H_

#ifdef __cplusplus
extern "C" {
#endif

// qore
#include "qore/QPSOLVER_DENSE/include/qpsolver_dense.h"
// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"

typedef struct dense_qp_qore_opts_
{
    int nsmax;       // maximum size of Schur complement
    int print_freq;  // print frequency,
                     // prtfreq  < 0: disable printing;
                     // prtfreq == 0: print on each call and include working set changes;
                     // prtfreq  > 0: print on every prtfreq seconds, but do not include working set
                     // changes;
    int warm_start;  // warm start with updated matrices H and C
    int warm_strategy;  // 0: ramp-up from zero homotopy; 1: setup homotopy from the previous
                        // solution
    int hot_start;      // hot start with unchanged matrices H and C
    int max_iter;       // maximum number of iterations
    int compute_t;      // compute t in qp_out (to have correct residuals in NLP)
} dense_qp_qore_opts;

typedef struct dense_qp_qore_memory_
{
    double *H;
    double *HH;
    double *g;
    double *gg;
    double *Zl;
    double *Zu;
    double *zl;
    double *zu;
    double *A;
    double *b;
    double *C;
    double *CC;
    double *Ct;
    double *CCt;
    double *d_lb0;
    double *d_ub0;
    double *d_lb;
    double *d_ub;
    double *d_lg;
    double *d_ug;
    double *d_ls;
    double *d_us;
    double *lb;
    double *ub;
    int *idxb;
    int *idxb_stacked;
    int *idxs;
    double *prim_sol;
    double *dual_sol;
    QoreProblemDense *QP;
    int num_iter;
    dense_qp_in *qp_stacked;
    double time_qp_solver_call;
    int iter;

} dense_qp_qore_memory;

acados_size_t dense_qp_qore_opts_calculate_size(void *config, dense_qp_dims *dims);
//
void *dense_qp_qore_opts_assign(void *config, dense_qp_dims *dims, void *raw_memory);
//
void dense_qp_qore_opts_initialize_default(void *config, dense_qp_dims *dims, void *opts_);
//
void dense_qp_qore_opts_update(void *config, dense_qp_dims *dims, void *opts_);
//
acados_size_t dense_qp_qore_memory_calculate_size(void *config, dense_qp_dims *dims, void *opts_);
//
void *dense_qp_qore_memory_assign(void *config, dense_qp_dims *dims, void *opts_, void *raw_memory);
//
acados_size_t dense_qp_qore_workspace_calculate_size(void *config, dense_qp_dims *dims, void *opts_);
//
int dense_qp_qore(void *config, dense_qp_in *qp_in, dense_qp_out *qp_out, void *opts_, void *memory_, void *work_);
//
void dense_qp_qore_eval_forw_sens(void *config_, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_);
//
void dense_qp_qore_memory_reset(void *config, void *qp_in, void *qp_out, void *opts, void *mem, void *work);
//
void dense_qp_qore_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2);
//
void dense_qp_qore_config_initialize_default(void *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_DENSE_QP_DENSE_QP_QORE_H_
