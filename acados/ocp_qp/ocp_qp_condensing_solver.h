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

#ifndef ACADOS_OCP_QP_OCP_QP_CONDENSING_SOLVER_H_
#define ACADOS_OCP_QP_OCP_QP_CONDENSING_SOLVER_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/ocp_qp/ocp_qp_condensing.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"


typedef struct ocp_qp_condensing_solver_args_ {
    ocp_qp_condensing_args *cond_args;
    dense_qp_solver *solver;
    void *solver_args;
} ocp_qp_condensing_solver_args;



typedef struct ocp_qp_condensing_solver_memory_ {
    ocp_qp_condensing_memory *cond_memory;
    void *solver_memory;
    dense_qp_in *qpd_in;
    dense_qp_out *qpd_out;
} ocp_qp_condensing_solver_memory;



typedef struct ocp_qp_condensing_solver_workspace_ {
    void *cond_work;
    void *solver_workspace;
    // TODO(dimitris): move from memory to workspace
    // dense_qp_in *qpd_in;
    // dense_qp_out *qpd_out;
} ocp_qp_condensing_solver_workspace;



//
int ocp_qp_condensing_solver_calculate_args_size(ocp_qp_dims *dims, void *solver_);
//
void *ocp_qp_condensing_solver_assign_args(ocp_qp_dims *dims, void *solver, void *raw_memory);
//
void ocp_qp_condensing_solver_initialize_default_args(void *args_);
//
int ocp_qp_condensing_solver_calculate_memory_size(ocp_qp_dims *dims, void *args_);
//
void *ocp_qp_condensing_solver_assign_memory(ocp_qp_dims *dims, void *args_, void *raw_memory);
//
int ocp_qp_condensing_solver_calculate_workspace_size(ocp_qp_dims *dims, void *args_);
//
int ocp_qp_condensing_solver(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *mem_, void *work_);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_CONDENSING_SOLVER_H_
