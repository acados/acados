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

#ifndef ACADOS_OCP_QP_OCP_QP_PARTIAL_CONDENSING_SOLVER_H_
#define ACADOS_OCP_QP_OCP_QP_PARTIAL_CONDENSING_SOLVER_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing.h"
#include "acados/utils/types.h"

typedef struct ocp_qp_partial_condensing_solver_opts_
{
    ocp_qp_partial_condensing_opts *pcond_opts;
    void *qp_solver_opts;
} ocp_qp_partial_condensing_solver_opts;

typedef struct ocp_qp_partial_condensing_solver_memory_
{
    ocp_qp_partial_condensing_memory *pcond_memory;
    void *solver_memory;
    ocp_qp_in *pcond_qp_in;
    ocp_qp_out *pcond_qp_out;
} ocp_qp_partial_condensing_solver_memory;

typedef struct ocp_qp_partial_condensing_solver_workspace_
{
    void *pcond_work;
    void *solver_work;
    // TODO(dimitris): move from memory to workspace
    // ocp_qp_in *pcond_qp_in;
    // ocp_qp_out *pcond_qp_out;
} ocp_qp_partial_condensing_solver_workspace;

//
int ocp_qp_partial_condensing_solver_opts_calculate_size(void *config, ocp_qp_dims *dims);
//
void *ocp_qp_partial_condensing_solver_opts_assign(void *config, ocp_qp_dims *dims,
                                                   void *raw_memory);
//
void ocp_qp_partial_condensing_solver_opts_initialize_default(void *config, ocp_qp_dims *dims,
                                                              void *opts_);
//
void ocp_qp_partial_condensing_solver_opts_update(void *config, ocp_qp_dims *dims, void *opts_);
//
void ocp_qp_partial_condensing_solver_opts_set(void *config_, void *opts_,
                                            const char *field, const void* value);
//
int ocp_qp_partial_condensing_solver_calculate_memory_size(void *config, ocp_qp_dims *dims,
                                                           void *opts_);
//
void *ocp_qp_partial_condensing_solver_memory_assign(void *config, ocp_qp_dims *dims, void *opts_,
                                                     void *raw_memory);
//
int ocp_qp_partial_condensing_solver_workspace_calculate_size(void *config, ocp_qp_dims *dims,
                                                              void *opts_);
//
int ocp_qp_partial_condensing_solver(void *config, ocp_qp_in *qp_in, ocp_qp_out *qp_out,
                                     void *opts_, void *mem_, void *work_);
//
void ocp_qp_partial_condensing_solver_config_initialize_default(void *config_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_PARTIAL_CONDENSING_SOLVER_H_
