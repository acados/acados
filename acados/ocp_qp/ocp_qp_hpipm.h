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

#ifndef ACADOS_OCP_QP_OCP_QP_HPIPM_H_
#define ACADOS_OCP_QP_OCP_QP_HPIPM_H_

#ifdef __cplusplus
extern "C" {
#endif

// hpipm
#include "hpipm/include/hpipm_d_ocp_qp_ipm.h"
// acados
#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"



// struct of arguments to the solver
typedef struct ocp_qp_hpipm_args_ {
    struct d_ocp_qp_ipm_arg *hpipm_args;
} ocp_qp_hpipm_args;



// struct of the solver memory
typedef struct ocp_qp_hpipm_memory_ {
    struct d_ocp_qp_ipm_workspace *hpipm_workspace;
} ocp_qp_hpipm_memory;

//
int ocp_qp_hpipm_calculate_args_size(ocp_qp_dims *dims);
//
void *ocp_qp_hpipm_assign_args(ocp_qp_dims *dims, void *raw_memory);
//
void ocp_qp_hpipm_initialize_default_args(void *args_);
//
int ocp_qp_hpipm_calculate_memory_size(ocp_qp_dims *dims, void *args_);
//
void *ocp_qp_hpipm_assign_memory(ocp_qp_dims *dims, void *args_, void *raw_memory);
//
int ocp_qp_hpipm_calculate_workspace_size(ocp_qp_dims *dims, void *args_);
//
int ocp_qp_hpipm(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *mem_, void *work_);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_HPIPM_H_
