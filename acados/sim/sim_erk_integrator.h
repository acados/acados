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
    int nx;
    int nu;
    int nz;
} sim_erk_dims;



typedef struct
{
    /* external functions */
    // explicit ode
    external_function_generic *expl_ode_fun;
    // hessian explicit ode
    external_function_generic *expl_ode_hes;
    // forward explicit vde
    external_function_generic *expl_vde_for;
    // adjoint explicit vde
    external_function_generic *expl_vde_adj;

} erk_model;



typedef struct
{
    // no memory
    void *dummy;
} sim_erk_memory;



typedef struct
{
    double *rhs_forw_in;  // x + S + p

    double *K_traj;         // (stages *nX) or (steps*stages*nX) for adj
    double *out_forw_traj;  // S or (steps+1)*nX for adj

    double *rhs_adj_in;
    double *out_adj_tmp;
    double *adj_traj;

} sim_erk_workspace;



// dims
int sim_erk_dims_calculate_size();
void *sim_erk_dims_assign(void *config_, void *raw_memory);
void sim_erk_dims_set(void *config_, void *dims_, const char *field, const int* value);
void sim_erk_dims_get(void *config_, void *dims_, const char *field, int* value);

// model
int sim_erk_model_calculate_size(void *config, void *dims);
void *sim_erk_model_assign(void *config, void *dims, void *raw_memory);
int sim_erk_model_set(void *model, const char *field, void *value);

// opts
int sim_erk_opts_calculate_size(void *config, void *dims);
//
void sim_erk_opts_update(void *config_, void *dims, void *opts_);
//
void *sim_erk_opts_assign(void *config, void *dims, void *raw_memory);
//
void sim_erk_opts_initialize_default(void *config, void *dims, void *opts_);
//
int sim_erk_opts_set(void *config_, void *opts_, const char *field, void *value);


// memory
int sim_erk_memory_calculate_size(void *config, void *dims, void *opts_);
//
void *sim_erk_memory_assign(void *config, void *dims, void *opts_, void *raw_memory);

// workspace
int sim_erk_workspace_calculate_size(void *config, void *dims, void *opts_);

//
int sim_erk(void *config, sim_in *in, sim_out *out, void *opts_, void *mem_, void *work_);
//
void sim_erk_config_initialize_default(void *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_ERK_INTEGRATOR_H_
