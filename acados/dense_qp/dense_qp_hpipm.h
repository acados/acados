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


#ifndef ACADOS_DENSE_QP_DENSE_QP_HPIPM_H_
#define ACADOS_DENSE_QP_DENSE_QP_HPIPM_H_

#ifdef __cplusplus
extern "C" {
#endif

// hpipm
#include "hpipm_d_dense_qp.h"
#include "hpipm_d_dense_qp_sol.h"
#include "hpipm_d_dense_qp_ipm.h"
// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"

typedef struct dense_qp_hpipm_args_ {
    struct d_dense_qp_ipm_arg *hpipm_args;
} dense_qp_hpipm_args;

typedef struct dense_qp_hpipm_memory_ {
    struct d_dense_qp_ipm_workspace *hpipm_workspace;
} dense_qp_hpipm_memory;

//
int dense_qp_hpipm_calculate_args_size(dense_qp_dims *dims);
//
void dense_qp_hpipm_initialize_default_args(void *args_);
//
void *dense_qp_hpipm_assign_args(dense_qp_dims *dims, void *raw_memory);
//
int dense_qp_hpipm_calculate_memory_size(dense_qp_dims *dims, void *args_);
//
void *dense_qp_hpipm_assign_memory(dense_qp_dims *dims, void *args_, void *raw_memory);
//
int dense_qp_hpipm_calculate_workspace_size(dense_qp_dims *dims, void *args_);
//
int dense_qp_hpipm(dense_qp_in *qp_in, dense_qp_out *qp_out, void *args_, void *mem_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_DENSE_QP_DENSE_QP_HPIPM_H_
