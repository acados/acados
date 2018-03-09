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


#ifndef ACADOS_C_LEGACY_CREATE_H_
#define ACADOS_C_LEGACY_CREATE_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include <acados/ocp_nlp/ocp_nlp_common.h>
#include <acados/ocp_nlp/ocp_nlp_gn_sqp.h>
#include <acados/ocp_qp/ocp_qp_full_condensing.h>
#include <acados/ocp_qp/ocp_qp_partial_condensing.h>
#include <acados/sim/sim_irk_integrator.h>
// acados_c
#include "acados_c/ocp_qp_interface.h"
#include "acados_c/sim_interface.h"

// ocp_qp_res *create_ocp_qp_res(ocp_qp_dims *dims);

// ocp_qp_res_ws *create_ocp_qp_res_ws(ocp_qp_dims *dims);

ocp_qp_full_condensing_args *ocp_qp_full_condensing_create_arguments(ocp_qp_dims *dims);

ocp_qp_full_condensing_memory *ocp_qp_full_condensing_create_memory(ocp_qp_dims *dims, ocp_qp_full_condensing_args *args);

ocp_qp_partial_condensing_args *ocp_qp_partial_condensing_create_arguments(ocp_qp_dims *dims);

ocp_qp_partial_condensing_memory *ocp_qp_partial_condensing_create_memory(ocp_qp_dims *dims, ocp_qp_partial_condensing_args *args);

ocp_nlp_in *create_ocp_nlp_in(ocp_nlp_solver_config *config, ocp_nlp_dims *dims);

ocp_nlp_out *create_ocp_nlp_out(ocp_nlp_solver_config *config, ocp_nlp_dims *dims);

ocp_nlp_gn_sqp_opts *ocp_nlp_gn_sqp_create_args(ocp_nlp_dims *dims, ocp_qp_solver_t qp_solver_name, sim_solver_t *sim_solver_names);

ocp_nlp_gn_sqp_memory *ocp_nlp_gn_sqp_create_memory(ocp_nlp_dims *dims, ocp_nlp_gn_sqp_opts *args);

sim_rk_opts *create_sim_irk_opts(sim_dims *dims);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_C_LEGACY_CREATE_H_
