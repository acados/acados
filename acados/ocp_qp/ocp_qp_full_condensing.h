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

#ifndef ACADOS_OCP_QP_OCP_QP_FULL_CONDENSING_H_
#define ACADOS_OCP_QP_OCP_QP_FULL_CONDENSING_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"



typedef struct ocp_qp_full_condensing_opts_
{
    struct d_cond_qp_ocp2dense_arg *hpipm_opts;
    int cond_hess; // 0 cond only rhs, 1 cond hess + rhs
    int expand_dual_sol; // 0 primal sol only, 1 primal + dual sol
	int ric_alg;
} ocp_qp_full_condensing_opts;



typedef struct ocp_qp_full_condensing_memory_
{
    struct d_cond_qp_ocp2dense_workspace *hpipm_workspace;
    // NOTE(dimitris): points to qp_in, does NOT copy to memory (needed for expansion)
    ocp_qp_in *qp_in;
} ocp_qp_full_condensing_memory;



//
void compute_dense_qp_dims(ocp_qp_dims *dims, dense_qp_dims *ddims);
//
int ocp_qp_full_condensing_opts_calculate_size(ocp_qp_dims *dims);
//
void *ocp_qp_full_condensing_opts_assign(ocp_qp_dims *dims, void *raw_memory);
//
void ocp_qp_full_condensing_opts_initialize_default(ocp_qp_dims *dims, void *opts_);
//
void ocp_qp_full_condensing_opts_update(ocp_qp_dims *dims, void *opts_);
//
void ocp_qp_full_condensing_opts_set(void *opts_, const char *field, void* value);
//
int ocp_qp_full_condensing_memory_calculate_size(ocp_qp_dims *dims, void *opts_);
//
void *ocp_qp_full_condensing_memory_assign(ocp_qp_dims *dims, void *opts_, void *raw_memory);
//
int ocp_qp_full_condensing_workspace_calculate_size(ocp_qp_dims *dims, void *opts_);
//
void ocp_qp_full_condensing(ocp_qp_in *in, dense_qp_in *out, ocp_qp_full_condensing_opts *opts,
                            ocp_qp_full_condensing_memory *mem, void *work);
//
void ocp_qp_full_expansion(dense_qp_out *in, ocp_qp_out *out, ocp_qp_full_condensing_opts *opts,
                           ocp_qp_full_condensing_memory *mem, void *work);
//
void ocp_qp_full_condensing_config_initialize_default(void *config_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_FULL_CONDENSING_H_
