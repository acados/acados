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
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_ipm.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"



typedef struct dense_qp_hpipm_opts_
{
    struct d_dense_qp_ipm_arg *hpipm_opts;
} dense_qp_hpipm_opts;



typedef struct dense_qp_hpipm_memory_
{
    struct d_dense_qp_ipm_workspace *hpipm_workspace;
} dense_qp_hpipm_memory;



//
int dense_qp_hpipm_opts_calculate_size(void *config, void *dims);
//
void *dense_qp_hpipm_opts_assign(void *config, void *dims, void *raw_memory);
//
void dense_qp_hpipm_opts_initialize_default(void *config, void *dims, void *opts_);
//
void dense_qp_hpipm_opts_update(void *config, void *dims, void *opts_);
//
int dense_qp_hpipm_calculate_memory_size(void *dims, void *opts_);
//
void *dense_qp_hpipm_assign_memory(void *dims, void *opts_, void *raw_memory);
//
int dense_qp_hpipm_calculate_workspace_size(void *dims, void *opts_);
//
int dense_qp_hpipm(void *config, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_);
//
void dense_qp_hpipm_config_initialize_default(void *config_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_DENSE_QP_DENSE_QP_HPIPM_H_
