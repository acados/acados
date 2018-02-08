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
#include "acados/utils/external_function.h"
#include "acados/utils/types.h"



typedef struct {
    external_function_fcn_ptrs *forward_vde;
    external_function_fcn_ptrs *adjoint_vde;
    external_function_fcn_ptrs *hess_vde;
} sim_erk_integrator_submodules;



typedef struct {
    // Options
    double interval;
    int num_stages;

    int num_steps;
    int num_forw_sens;

    double *A_mat;
    double *c_vec;
    double *b_vec;

    bool sens_forw;
    bool sens_adj;
    bool sens_hess;

    // Submodules
    sim_erk_integrator_submodules submodules;

    // Arguments for functions
    void *forward_vde_args;
    void *adjoint_vde_args;
    void *hess_vde_args;
} sim_erk_integrator_args;



typedef struct {
    // Memory for functions
    void *forward_vde_mem;
    void *adjoint_vde_mem;
    void *hess_vde_mem;
} sim_erk_integrator_memory;



typedef struct {
    double *rhs_forw_in;  // x + S + p

    double *K_traj; // (stages *nX) or (steps*stages*nX) for adj
    double *out_forw_traj; // S or (steps+1)*nX for adj

    double *rhs_adj_in;
    double *out_adj_tmp;
    double *adj_traj;

    // Workspace for functions
    void *forward_vde_work;
    void *adjoint_vde_work;
    void *hess_vde_work;
} sim_erk_integrator_workspace;



//
int sim_erk_integrator_calculate_args_size(sim_dims *dims, void *submodules_);
//
void *sim_erk_integrator_assign_args(sim_dims *dims, void **submodules_, void *raw_memory);
//
void *sim_erk_integrator_copy_args(sim_dims *dims, void *raw_memory, void *source_);
//
void sim_erk_integrator_initialize_default_args(sim_dims *dims, void *args_);
//
int sim_erk_integrator_calculate_memory_size(sim_dims *dims, void *args_);
//
void *sim_erk_integrator_assign_memory(sim_dims *dims, void *args_, void *raw_memory);
//
int sim_erk_integrator_calculate_workspace_size(sim_dims *dims, void *args_);
//
int sim_erk_integrator(sim_in *in, sim_out *out, void *args_, void *mem_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_ERK_INTEGRATOR_H_
