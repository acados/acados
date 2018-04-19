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

#ifndef INTERFACES_ACADOS_C_DENSE_QP_INTERFACE_H_
#define INTERFACES_ACADOS_C_DENSE_QP_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/dense_qp/dense_qp_common.h"

typedef enum { DENSE_QP_HPIPM, DENSE_QP_QORE, DENSE_QP_QPOASES } dense_qp_solver_t;

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
int dense_qp_calculate_size(qp_solver_config *config, dense_qp_dims *dims, void *opts_);
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

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // INTERFACES_ACADOS_C_DENSE_QP_INTERFACE_H_
