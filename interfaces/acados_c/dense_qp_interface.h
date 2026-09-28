/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef INTERFACES_ACADOS_C_DENSE_QP_INTERFACE_H_
#define INTERFACES_ACADOS_C_DENSE_QP_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/dense_qp/dense_qp_common.h"

typedef enum { DENSE_QP_HPIPM, DENSE_QP_QORE, DENSE_QP_QPOASES, DENSE_QP_OOQP, DENSE_QP_DAQP } dense_qp_solver_t;

typedef struct
{
    dense_qp_solver_t qp_solver;
} dense_qp_solver_plan;

typedef struct
{
    qp_solver_config *config;
    void *dims;
    void *opts;
    void *mem;
    void *work;
} dense_qp_solver;

qp_solver_config *dense_qp_config_create(dense_qp_solver_plan *plan);
//
dense_qp_dims *dense_qp_dims_create();
//
dense_qp_in *dense_qp_in_create(qp_solver_config *config, dense_qp_dims *dims);
//
dense_qp_out *dense_qp_out_create(qp_solver_config *config, dense_qp_dims *dims);
//
void *dense_qp_opts_create(qp_solver_config *config, dense_qp_dims *dims);
//
acados_size_t dense_qp_calculate_size(qp_solver_config *config, dense_qp_dims *dims, void *opts_);
//
dense_qp_solver *dense_qp_assign(qp_solver_config *config, dense_qp_dims *dims, void *opts_,
                                 void *raw_memory);
//
dense_qp_solver *dense_qp_create(qp_solver_config *config, dense_qp_dims *dims, void *opts_);
//
int dense_qp_solve(dense_qp_solver *solver, dense_qp_in *qp_in, dense_qp_out *qp_out);
//
void dense_qp_inf_norm_residuals(dense_qp_dims *dims, dense_qp_in *qp_in, dense_qp_out *qp_out,
                                 double *res);
//
bool dense_qp_set_field_double_array(const char *field, double *arr, dense_qp_in *qp_in);
//
bool dense_qp_set_field_int_array(const char *field, int *arr, dense_qp_in *qp_in);
//
bool dense_qp_get_field_double_array(const char *field, dense_qp_in *qp_in, double *arr);
//
bool dense_qp_get_field_int_array(const char *field, dense_qp_in *qp_in, int *arr);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // INTERFACES_ACADOS_C_DENSE_QP_INTERFACE_H_
