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

#ifndef ACADOS_UTILS_CASADI_WRAPPER_H_
#define ACADOS_UTILS_CASADI_WRAPPER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"

typedef struct {
    int dummy;
} casadi_wrapper_dims;

typedef struct {
    double *x;
    double *u;
    double *p;
    double *mul;

    bool compute_y;
    bool compute_jac_y;
    bool compute_grad_mul_y;
    bool compute_hess_mul_y;
} casadi_wrapper_in;

typedef struct {
    double *y;
    double *jac_y;
    double *grad_mul_y;
    double *hess_mul_y;
} casadi_wrapper_out;

typedef struct {
    int_t (*fun)(const real_t **, real_t **, int_t *, real_t *, int_t);
    int_t (*dims)(int_t *, int_t *, int_t *, int_t *);
    const int_t *(*sparsity)(int_t);
} casadi_wrapper_args;

typedef struct {
    int dummy;
} casadi_wrapper_memory;

typedef struct {
    const real_t **arg;
    real_t **res;
    int_t *iw;
    real_t *w;
    real_t **sparse_res;
} casadi_wrapper_workspace;

//
int casadi_wrapper_calculate_args_size(casadi_wrapper_dims *dims);
//
void *casadi_wrapper_assign_args(casadi_wrapper_dims *dims, void *raw_memory);
//
void casadi_wrapper_initialize_default_args(casadi_wrapper_args *args);
//
int casadi_wrapper_calculate_memory_size(casadi_wrapper_dims *dims, casadi_wrapper_args *args);
//
void *casadi_wrapper_assign_memory(casadi_wrapper_dims *dims, casadi_wrapper_args *args, void *raw_memory);
//
int casadi_wrapper_calculate_workspace_size(casadi_wrapper_dims *dims, casadi_wrapper_args *args);
//
int casadi_wrapper(casadi_wrapper_in *cw_in, casadi_wrapper_out *cw_out, casadi_wrapper_args *args, casadi_wrapper_memory *mem, casadi_wrapper_workspace *work);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_CASADI_WRAPPER_H_
