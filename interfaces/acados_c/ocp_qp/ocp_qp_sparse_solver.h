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

#ifndef ACADOS_C_OCP_QP_OCP_QP_SPARSE_SOLVER_H_
#define ACADOS_C_OCP_QP_OCP_QP_SPARSE_SOLVER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <acados/ocp_qp/ocp_qp_sparse_solver.h>

#include "acados_c/ocp_qp.h"

//
// void *ocp_qp_sparse_solver_copy_args(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *raw_memory, void *source_);
//
int ocp_qp_sparse_solver_calculate_submodules_size(ocp_qp_solver_config *config, ocp_qp_dims *dims);
//
void *ocp_qp_sparse_solver_assign_submodules(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *raw_memory);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_C_OCP_QP_OCP_QP_SPARSE_SOLVER_H_