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

#define HPIPM_SP


// hpipm
#include "hpipm/include/hpipm_d_ocp_qp_ipm.h"

#ifdef HPIPM_SP
#include "hpipm/include/hpipm_s_ocp_qp_ipm.h"
#endif

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

// struct of arguments to the solver
// TODO(roversch): why not make this a typedef of the underlying struct?
#ifdef HPIPM_SP
typedef struct ocp_qp_hpipm_opts_
{
    struct s_ocp_qp_ipm_arg *hpipm_opts;
} ocp_qp_hpipm_opts;
#else
typedef struct ocp_qp_hpipm_opts_
{
    struct d_ocp_qp_ipm_arg *hpipm_opts;
} ocp_qp_hpipm_opts;
#endif

// TODO(roversch): why not make this a typedef of the underlying struct?
// struct of the solver memory
#ifdef HPIPM_SP
typedef struct ocp_qp_hpipm_memory_
{
    struct s_ocp_qp_ipm_workspace *hpipm_workspace;
} ocp_qp_hpipm_memory;
#else
typedef struct ocp_qp_hpipm_memory_
{
    struct d_ocp_qp_ipm_workspace *hpipm_workspace;
} ocp_qp_hpipm_memory;
#endif
//
int ocp_qp_hpipm_opts_calculate_size(void *config, void *dims);
//
void *ocp_qp_hpipm_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_qp_hpipm_opts_initialize_default(void *config, void *dims, void *opts_);
//
void ocp_qp_hpipm_opts_update(void *config, void *dims, void *opts_);
//
int ocp_qp_hpipm_memory_calculate_size(void *config, void *dims, void *opts_);
//
void *ocp_qp_hpipm_memory_assign(void *config, void *dims, void *opts_, void *raw_memory);
//
int ocp_qp_hpipm_workspace_calculate_size(void *config, void *dims, void *opts_);
//
int ocp_qp_hpipm(void *config, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_);
//
void ocp_qp_hpipm_config_initialize_default(void *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_HPIPM_H_
