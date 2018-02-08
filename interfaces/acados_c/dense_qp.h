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

#ifndef ACADOS_C_DENSE_QP_H_
#define ACADOS_C_DENSE_QP_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include <acados/dense_qp/dense_qp_common.h>
// acados_c
#include "acados_c/common.h"

typedef enum {
    DENSE_QP_HPIPM,
    DENSE_QP_QORE,
    DENSE_QP_QPOASES
} dense_qp_solver_t;

typedef struct {
    dense_qp_solver_t qp_solver;
} dense_qp_solver_config;

typedef struct {
    dense_qp_solver_fcn_ptrs *fcn_ptrs;
    void *dims;
    void *args;
    void *mem;
    void *work;
} dense_qp_solver;

// INPUT, OUTPUT AND OPTIONS
//
dense_qp_dims *create_dense_qp_dims();
//
dense_qp_in *create_dense_qp_in(dense_qp_dims *dims);
//
dense_qp_out *create_dense_qp_out(dense_qp_dims *dims);
//
int dense_qp_calculate_args_size(dense_qp_solver_fcn_ptrs *fcn_ptrs, dense_qp_dims *dims);
//
void *dense_qp_assign_args(dense_qp_solver_fcn_ptrs *fcn_ptrs, dense_qp_dims *dims, void *raw_memory);
//
void *dense_qp_create_args(dense_qp_solver_fcn_ptrs *fcn_ptrs, dense_qp_dims *dims);
//
void *dense_qp_copy_args(dense_qp_solver_fcn_ptrs *fcn_ptrs, dense_qp_dims *dims, void *raw_memory, void *source);

// BASIC INTERFACE
//
int dense_qp_calculate_size(dense_qp_solver_fcn_ptrs *fcn_ptrs, dense_qp_dims *dims, void *args_);
//
dense_qp_solver *dense_qp_assign(dense_qp_solver_fcn_ptrs *fcn_ptrs, dense_qp_dims *dims, void *args_, void *raw_memory);
//
dense_qp_solver *dense_qp_create(dense_qp_solver_fcn_ptrs *fcn_ptrs, dense_qp_dims *dims, void *args_);
//
int dense_qp_solve(dense_qp_solver *solver, dense_qp_in *qp_in, dense_qp_out *qp_out);

// OPTIONS BASED CONFIGURATION STRATEGY
//
int dense_qp_calculate_submodules_size(dense_qp_solver_config *config, dense_qp_dims *dims);
//
void *dense_qp_assign_submodules(dense_qp_solver_config *config, dense_qp_dims *dims, void *raw_memory);
//
int calculate_dense_qp_solver_fcn_ptrs_size(dense_qp_solver_config *config, dense_qp_dims *dims);
//
void *assign_dense_qp_solver_fcn_ptrs(dense_qp_solver_config *config, dense_qp_dims *dims, void *raw_memory);
//
void *create_dense_qp_solver_fcn_ptrs(dense_qp_solver_config *config, dense_qp_dims *dims);

// EXPERT INTERFACE
//
int set_dense_qp_solver_fcn_ptrs(dense_qp_solver_config *config, dense_qp_solver_fcn_ptrs *fcn_ptrs);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_C_DENSE_QP_H_
