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
#include "acados/utils/external_function.h"



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
int casadi_wrapper_calculate_args_size(external_function_dims *dims, void *submodules_);
//
void *casadi_wrapper_assign_args(external_function_dims *dims, void **submodules_, void *raw_memory);
//
void *casadi_wrapper_copy_args(external_function_dims *dims, void *raw_memory, void *source_);
//
void casadi_wrapper_initialize_default_args(void *args_);
//
int casadi_wrapper_calculate_memory_size(external_function_dims *dims, void *args_);
//
void *casadi_wrapper_assign_memory(external_function_dims *dims, void *args_, void *raw_memory);
//
int casadi_wrapper_calculate_workspace_size(external_function_dims *dims, void *args_);
//
int casadi_wrapper(external_function_in *ef_in, external_function_out *ef_out, void *args_, void *mem_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_CASADI_WRAPPER_H_
