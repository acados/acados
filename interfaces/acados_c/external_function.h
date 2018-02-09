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

#ifndef ACADOS_C_EXTERNAL_FUNCTION_H_
#define ACADOS_C_EXTERNAL_FUNCTION_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include <acados/utils/external_function.h>
// acados_c
#include "acados_c/common.h"

typedef enum {
    CASADI_WRAPPER
} external_function_t;

typedef struct {
    external_function_t type;
} external_function_config;

typedef struct {
    external_function_fcn_ptrs *fcn_ptrs;
    void *dims;
    void *args;
    void *mem;
    void *work;
} external_function;

// INPUT, OUTPUT AND OPTIONS
//
external_function_dims *create_external_function_dims(int num_inputs, int num_outputs);
//
external_function_in *create_external_function_in(external_function_dims *dims);
//
external_function_out *create_external_function_out(external_function_dims *dims);
//
int external_function_calculate_args_size(external_function_fcn_ptrs *fcn_ptrs, external_function_dims *dims);
//
void *external_function_assign_args(external_function_fcn_ptrs *fcn_ptrs, external_function_dims *dims, void *raw_memory);
//
void *external_function_create_args(external_function_fcn_ptrs *fcn_ptrs, external_function_dims *dims);
//
void *external_function_copy_args(external_function_fcn_ptrs *fcn_ptrs, external_function_dims *dims, void *raw_memory, void *source);

// BASIC INTERFACE
//
int external_function_calculate_size(external_function_fcn_ptrs *fcn_ptrs, external_function_dims *dims, void *args_);
//
external_function *external_function_assign(external_function_fcn_ptrs *fcn_ptrs, external_function_dims *dims, void *args_, void *raw_memory);
//
external_function *external_function_create(external_function_fcn_ptrs *fcn_ptrs, external_function_dims *dims, void *args_);
//
int external_function_eval(external_function *ext_fun, external_function_in *ef_in, external_function_out *ef_out);

// OPTIONS BASED CONFIGURATION STRATEGY
//
int external_function_calculate_submodules_size(external_function_config *config, external_function_dims *dims);
//
void *external_function_assign_submodules(external_function_config *config, external_function_dims *dims, void *raw_memory);
//
int calculate_external_function_fcn_ptrs_size(external_function_config *config, external_function_dims *dims);
//
void *assign_external_function_fcn_ptrs(external_function_config *config, external_function_dims *dims, void *raw_memory);
//
void *create_external_function_fcn_ptrs(external_function_config *config, external_function_dims *dims);

// EXPERT INTERFACE
//
int set_external_function_fcn_ptrs(external_function_config *config, external_function_fcn_ptrs *fcn_ptrs);




#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_C_EXTERNAL_FUNCTION_H_