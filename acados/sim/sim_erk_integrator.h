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

#include "acados/sim/sim_common.h"
#include "acados/utils/types.h"

typedef struct {

} sim_erk_memory;

typedef struct {

    double *rhs_forw_in;  // x + S + p

    double *K_traj; // (stages *nX) or (steps*stages*nX) for adj
    double *out_forw_traj; // S or (steps+1)*nX for adj

    double *rhs_adj_in;
    double *out_adj_tmp;
    double *adj_traj;

} sim_erk_workspace;

int sim_erk_opts_calculate_size(sim_dims *dims);

void *sim_erk_assign_opts(sim_dims *dims, void *raw_memory);

void sim_erk_initialize_default_args(sim_dims *dims, void *opts_);

int sim_erk_calculate_memory_size(sim_dims *dims, void *opts_);

void *sim_erk_assign_memory(sim_dims *dims, void *opts_, void *raw_memory);

void *sim_erk_create_memory(sim_dims *dims, void *opts_);

int sim_erk(sim_in *in, sim_out *out, void *opts_, void *mem_, void *work_);

int sim_erk_calculate_workspace_size(sim_dims *dims, void *opts_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_ERK_INTEGRATOR_H_
