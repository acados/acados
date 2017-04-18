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

#ifndef ACADOS_SIM_SIM_ERK_INTEGRATOR_H_
#define ACADOS_SIM_SIM_ERK_INTEGRATOR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/sim/sim_rk_common.h"
#include "acados/utils/types.h"

typedef struct {
    real_t *K_traj;

    real_t *rhs_forw_in;
    real_t *out_forw_traj;

    real_t *adj_traj;
    real_t *rhs_adj_in;
    real_t *out_adj_tmp;
} sim_erk_workspace;


int_t sim_erk(const sim_in *in, sim_out *out, void *args, void *mem, void *work);

void sim_erk_create_workspace(const sim_in *in, void *args, sim_erk_workspace *work);

void sim_erk_create_arguments(void *args, int_t num_stages);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_ERK_INTEGRATOR_H_
