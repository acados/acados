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

#ifndef ACADOS_C_OCP_lin_H_
#define ACADOS_C_OCP_lin_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include <acados/ocp_lin/ocp_lin_common.h>
// acados_c
#include "acados_c/common.h"

typedef enum {
    GAUSS_NEWTON,
    EXACT_HESSIAN,
    DOPUS
} ocp_lin_method_t;

typedef struct {
    ocp_lin_method_t lin_method;
} ocp_lin_method_config;

typedef struct {
    ocp_lin_method_fcn_ptrs *fcn_ptrs;
    void *dims;
    void *args;
    void *mem;
    void *work;
} ocp_lin_method;

// INPUT, OUTPUT AND OPTIONS
//
ocp_lin_dims *create_ocp_lin_dims();
//
ocp_lin_in *create_ocp_lin_in(ocp_lin_dims *dims);
//
ocp_lin_out *create_ocp_lin_out(ocp_lin_dims *dims);
//
int ocp_lin_calculate_args_size(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims);
//
void *ocp_lin_assign_args(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims, void *raw_memory);
//
void *ocp_lin_create_args(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims);
//
void *ocp_lin_copy_args(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims, void *raw_memory, void *source);

// BASIC INTERFACE
//
int ocp_lin_calculate_size(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims, void *args_);
//
ocp_lin_method *ocp_lin_assign(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims, void *args_, void *raw_memory);
//
ocp_lin_method *ocp_lin_create(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims, void *args_);
//
int ocp_lin_solve(ocp_lin_method *solver, ocp_lin_in *qp_in, ocp_lin_out *qp_out);

// OPTIONS BASED CONFIGURATION STRATEGY
//
int ocp_lin_calculate_submodules_size(ocp_lin_method_config *config, ocp_lin_dims *dims);
//
void *ocp_lin_assign_submodules(ocp_lin_method_config *config, ocp_lin_dims *dims, void *raw_memory);
//
int calculate_ocp_lin_solver_fcn_ptrs_size(ocp_lin_method_config *config, ocp_lin_dims *dims);
//
void *assign_ocp_lin_solver_fcn_ptrs(ocp_lin_method_config *config, ocp_lin_dims *dims, void *raw_memory);
//
void *create_ocp_lin_solver_fcn_ptrs(ocp_lin_method_config *config, ocp_lin_dims *dims);

// EXPERT INTERFACE
//
int set_ocp_lin_method_fcn_ptrs(ocp_lin_method_config *config, ocp_lin_method_fcn_ptrs *fcn_ptrs);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_C_OCP_lin_QP_H_
