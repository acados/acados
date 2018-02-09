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

#ifndef ACADOS_UTILS_EXTERNAL_FUNCTION_H_
#define ACADOS_UTILS_EXTERNAL_FUNCTION_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"



typedef struct {
    int num_inputs;
    int num_outputs;
    int *input_dims;
    int *output_dims;
} external_function_dims;



typedef struct {
    double **inputs;
    bool *compute_output;
} external_function_in;



typedef struct {
    double **outputs;
} external_function_out;



typedef struct {
    int (*fun)(external_function_in *in, external_function_out *out, void *args, void *mem, void *work);
    int (*calculate_args_size)(external_function_dims *dims, void *submodules);
    void *(*assign_args)(external_function_dims *dims, void **submodules, void *raw_memory);
    void *(*copy_args)(external_function_dims *dims, void *raw_memory, void *source);
    void (*initialize_default_args)(void *args);
    int (*calculate_memory_size)(external_function_dims *dims, void *args);
    void *(*assign_memory)(external_function_dims *dims, void *args, void *raw_memory);
    int (*calculate_workspace_size)(external_function_dims *dims, void *args);
    void *submodules;
} external_function_fcn_ptrs;



//
int external_function_dims_calculate_size(int num_inputs, int num_outputs);
//
external_function_dims *assign_external_function_dims(int num_inputs, int num_outputs, void *raw_memory);
//
int external_function_in_calculate_size(external_function_dims *dims);
//
external_function_in *assign_external_function_in(external_function_dims *dims, void *raw_memory);
//
int external_function_out_calculate_size(external_function_dims *dims);
//
external_function_out *assign_external_function_out(external_function_dims *dims, void *raw_memory);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_EXTERNAL_FUNCTION_H_