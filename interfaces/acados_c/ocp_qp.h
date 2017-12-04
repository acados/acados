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
#include <acados/utils/types.h>

typedef enum {
    OCP_QP_HPIPM,
    OCP_QP_HPMPC,
    OCP_QP_OOQP,
    OCP_QP_QPDUNES,
    PARTIAL_CONDENSING_HPIPM,
    PARTIAL_CONDENSING_HPMPC,
    PARTIAL_CONDENSING_OOQP,
    PARTIAL_CONDENSING_QPDUNES,
    CONDENSING_HPIPM,
    CONDENSING_QPOASES,
    CONDENSING_QORE
} ocp_qp_solver_t;

typedef struct {
    ocp_qp_solver_t qp_solver;
} ocp_qp_config;

//
int ocp_qp_calculate_args_size(ocp_qp_config *config, ocp_qp_dims *dims);
//
void *ocp_qp_assign_args(ocp_qp_config *config, ocp_qp_dims *dims, void *raw_memory);
//
void *ocp_qp_create_args(ocp_qp_config *config, ocp_qp_dims *dims);
//
void ocp_qp_assign_default_args(ocp_qp_config *config, void *args_);
//
int ocp_qp_calculate_memory_size(ocp_qp_dims *dims, void *args_);
//
void *ocp_qp_assign_memory(ocp_qp_dims *dims, void *args_, void *raw_memory);
//
void *ocp_qp_create_memory(ocp_qp_dims *dims, void *args_);
//
int ocp_qp_calculate_workspace_size(ocp_qp_dims *dims, void *args_);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_C_OCP_QP_H_
