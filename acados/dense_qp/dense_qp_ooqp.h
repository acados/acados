/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef ACADOS_DENSE_QP_DENSE_QP_OOQP_H_
#define ACADOS_DENSE_QP_DENSE_QP_OOQP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"

enum dense_qp_ooqp_termination_code
{
  DENSE_SUCCESSFUL_TERMINATION = 0,
  DENSE_NOT_FINISHED,
  DENSE_MAX_ITS_EXCEEDED,
  DENSE_INFEASIBLE,
  DENSE_UNKNOWN
};

typedef struct dense_qp_ooqp_opts_
{
    int printLevel;
    int useDiagonalWeights;  // TODO(dimitris): implement option
    int fixHessian;
    int fixDynamics;
    int fixInequalities;
} dense_qp_ooqp_opts;

typedef struct dense_qp_ooqp_workspace_
{
    double *x;
    double *gamma;
    double *phi;
    double *y;
    double *z;
    double *lambda;
    double *pi;
    double objectiveValue;
} dense_qp_ooqp_workspace;

typedef struct dense_qp_ooqp_memory_
{
    int firstRun;
    int nx;
    int my;
    int mz;
    double *c;
    double *dQ;
    double *xlow;
    char *ixlow;
    double *xupp;
    char *ixupp;
    double *dA;
    double *bA;
    double *dC;
    double *clow;
    char *iclow;
    double *cupp;
    char *icupp;
    double time_qp_solver_call;
    int iter;

} dense_qp_ooqp_memory;

//
acados_size_t dense_qp_ooqp_opts_calculate_size(void *config_, dense_qp_dims *dims);
//
void *dense_qp_ooqp_opts_assign(void *config_, dense_qp_dims *dims, void *raw_memory);
//
void dense_qp_ooqp_opts_initialize_default(void *config_, dense_qp_dims *dims, void *opts_);
//
void dense_qp_ooqp_opts_update(void *config_, dense_qp_dims *dims, void *opts_);
//
acados_size_t dense_qp_ooqp_memory_calculate_size(void *config_, dense_qp_dims *dims, void *opts_);
//
void *dense_qp_ooqp_memory_assign(void *config_, dense_qp_dims *dims, void *opts_,
                                  void *raw_memory);
//
acados_size_t dense_qp_ooqp_workspace_calculate_size(void *config_, dense_qp_dims *dims, void *opts_);
//
int dense_qp_ooqp(void *config_, dense_qp_in *qp_in, dense_qp_out *qp_out, void *opts_,
                  void *memory_, void *work_);
//
void dense_qp_ooqp_destroy(void *mem_, void *work);
//
void dense_qp_ooqp_memory_reset(void *config, void *qp_in, void *qp_out, void *opts, void *mem, void *work);
//
void dense_qp_ooqp_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2);
//
void dense_qp_ooqp_config_initialize_default(void *config_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_DENSE_QP_DENSE_QP_OOQP_H_
