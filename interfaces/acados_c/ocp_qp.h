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

#ifndef ACADOS_C_OCP_QP_H_
#define ACADOS_C_OCP_QP_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include <acados/ocp_qp/ocp_qp_common.h>
#include <acados/ocp_qp/ocp_qp_common_frontend.h>
// acados_c
#include "acados_c/common.h"



typedef enum {
    SPARSE_QP_HPIPM,
    SPARSE_QP_HPMPC,
    SPARSE_QP_OOQP,
    SPARSE_QP_QPDUNES,
    PARTIAL_CONDENSING_HPIPM,
    PARTIAL_CONDENSING_HPMPC,
    PARTIAL_CONDENSING_OOQP,
    PARTIAL_CONDENSING_QPDUNES,
    FULL_CONDENSING_HPIPM,
    FULL_CONDENSING_QPOASES,
    FULL_CONDENSING_QORE
} ocp_qp_solver_t;



typedef struct {
    ocp_qp_solver_t qp_solver;
} ocp_qp_solver_config;



typedef struct {
    ocp_qp_solver_fcn_ptrs *fcn_ptrs;
    void *dims;
    void *args;
    void *mem;
    void *work;
} ocp_qp_solver;



// INPUT, OUTPUT AND OPTIONS
//
ocp_qp_dims *create_ocp_qp_dims(int N);
//
ocp_qp_in *create_ocp_qp_in(ocp_qp_dims *dims);
//
ocp_qp_out *create_ocp_qp_out(ocp_qp_dims *dims);
//
int ocp_qp_calculate_args_size(ocp_qp_solver_fcn_ptrs *fcn_ptrs, ocp_qp_dims *dims);
//
void *ocp_qp_assign_args(ocp_qp_solver_fcn_ptrs *fcn_ptrs, ocp_qp_dims *dims, void *raw_memory);
//
void *ocp_qp_create_args(ocp_qp_solver_fcn_ptrs *fcn_ptrs, ocp_qp_dims *dims);
//
void *ocp_qp_copy_args(ocp_qp_solver_fcn_ptrs *fcn_ptrs, ocp_qp_dims *dims, void *raw_memory, void *source);

// BASIC INTERFACE
//
int ocp_qp_calculate_size(ocp_qp_solver_fcn_ptrs *fcn_ptrs, ocp_qp_dims *dims, void *args_);
//
ocp_qp_solver *ocp_qp_assign(ocp_qp_solver_fcn_ptrs *fcn_ptrs, ocp_qp_dims *dims, void *args_, void *raw_memory);
//
ocp_qp_solver *ocp_qp_create(ocp_qp_solver_fcn_ptrs *fcn_ptrs, ocp_qp_dims *dims, void *args_);
//
int ocp_qp_solve(ocp_qp_solver *solver, ocp_qp_in *qp_in, ocp_qp_out *qp_out);

// OPTIONS BASED CONFIGURATION STRATEGY
//
int ocp_qp_calculate_submodules_size(ocp_qp_solver_config *config, ocp_qp_dims *dims);
//
void *ocp_qp_assign_submodules(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *raw_memory);
//
int calculate_ocp_qp_solver_fcn_ptrs_size(ocp_qp_solver_config *config, ocp_qp_dims *dims);
//
void *assign_ocp_qp_solver_fcn_ptrs(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *raw_memory);
//
void *create_ocp_qp_solver_fcn_ptrs(ocp_qp_solver_config *config, ocp_qp_dims *dims);

// EXPERT INTERFACE
//
int set_ocp_qp_solver_fcn_ptrs(ocp_qp_solver_config *config, ocp_qp_solver_fcn_ptrs *fcn_ptrs);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_C_OCP_QP_H_
