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

#ifndef ACADOS_OCP_QP_OCP_QP_PARTIAL_CONDENSING_H_
#define ACADOS_OCP_QP_OCP_QP_PARTIAL_CONDENSING_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/ocp_qp/ocp_qp_common.h"



typedef struct ocp_qp_partial_condensing_opts_
{
    struct d_cond_qp_ocp2ocp_arg *hpipm_opts;
    ocp_qp_dims *pcond_dims;  // TODO(all): move to dims
    int *block_size;
    int N2;
    int N2_bkp;
	int ric_alg;
} ocp_qp_partial_condensing_opts;



typedef struct ocp_qp_partial_condensing_memory_
{
    struct d_cond_qp_ocp2ocp_workspace *hpipm_workspace;
    ocp_qp_in *qp_in;
    ocp_qp_in *pcond_qp_in;
} ocp_qp_partial_condensing_memory;



//
int ocp_qp_partial_condensing_opts_calculate_size(ocp_qp_dims *dims);
//
void *ocp_qp_partial_condensing_opts_assign(ocp_qp_dims *dims, void *raw_memory);
//
void ocp_qp_partial_condensing_opts_initialize_default(ocp_qp_dims *dims, void *opts_);
//
void ocp_qp_partial_condensing_opts_update(ocp_qp_dims *dims, void *opts_);
//
void ocp_qp_partial_condensing_opts_set(void *opts_, const char *field, void* value);
//
int ocp_qp_partial_condensing_memory_calculate_size(ocp_qp_dims *dims, void *opts_);
//
void *ocp_qp_partial_condensing_memory_assign(ocp_qp_dims *dims, void *opts, void *raw_memory);
//
int ocp_qp_partial_condensing_workspace_calculate_size(ocp_qp_dims *dims, void *opts_);
//
void ocp_qp_partial_condensing(ocp_qp_in *in, ocp_qp_in *out, ocp_qp_partial_condensing_opts *opts,
                               ocp_qp_partial_condensing_memory *mem, void *work);
//
void ocp_qp_partial_expansion(ocp_qp_out *in, ocp_qp_out *out, ocp_qp_partial_condensing_opts *opts,
                              ocp_qp_partial_condensing_memory *mem, void *work);
//
void ocp_qp_partial_condensing_config_initialize_default(void *config_);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_PARTIAL_CONDENSING_H_
