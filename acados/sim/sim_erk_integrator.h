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



typedef struct
{
	/* external functions */
	// explicit ode
	external_function_generic *ode_expl;
	// jacobian explicit ode
	external_function_generic *jac_ode_expl;
	// hessian explicit ode
	external_function_generic *hess_ode_expl;
	// forward explicit vde
	external_function_generic *forw_vde_expl;
	// adjoint explicit vde
	external_function_generic *adj_vde_expl;

} erk_model;



typedef struct
{
	// no memory
} sim_erk_memory;



typedef struct
{

    double *rhs_forw_in;  // x + S + p

    double *K_traj; // (stages *nX) or (steps*stages*nX) for adj
    double *out_forw_traj; // S or (steps+1)*nX for adj

    double *rhs_adj_in;
    double *out_adj_tmp;
    double *adj_traj;

} sim_erk_workspace;



//
int sim_erk_model_calculate_size(void *config, sim_dims *dims);
//
void *sim_erk_model_assign(void *config, sim_dims *dims, void *raw_memory);
//
void sim_erk_model_set_forward_vde(sim_in *in, void *fun);
//
void sim_erk_model_set_adjoint_vde(sim_in *in, void *fun);

int sim_erk_opts_calculate_size(void *config, sim_dims *dims);

void *sim_erk_opts_assign(void *config, sim_dims *dims, void *raw_memory);

void sim_erk_opts_initialize_default(void *config, sim_dims *dims, void *opts_);

int sim_erk_memory_calculate_size(void *config, sim_dims *dims, void *opts_);

void *sim_erk_memory_assign(void *config, sim_dims *dims, void *opts_, void *raw_memory);

int sim_erk_workspace_calculate_size(void *config, sim_dims *dims, void *opts_);

int sim_erk(void *config, sim_in *in, sim_out *out, void *opts_, void *mem_, void *work_);

void sim_erk_config_initialize_default(void *config);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_ERK_INTEGRATOR_H_
