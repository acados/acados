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

#ifndef ACADOS_SIM_SIM_LIFTED_IRK_INTEGRATOR_H_
#define ACADOS_SIM_SIM_LIFTED_IRK_INTEGRATOR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/sim/sim_collocation.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/external_function.h"
#include "acados/utils/types.h"

#define TRIPLE_LOOP 1
#define CODE_GENERATION 0



typedef struct {
    external_function_fcn_ptrs *forward_vde;
    external_function_fcn_ptrs *jacobian_ode;
} sim_lifted_irk_integrator_submodules;



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

    int newton_iter;
    Newton_scheme *scheme;

    // Submodules
    sim_lifted_irk_integrator_submodules submodules;

    // Arguments for functions
    void *forward_vde_args;
    void *jacobian_ode_args;
} sim_lifted_irk_integrator_args;



typedef struct {

    double *grad_correction;
    double *grad_K;  // gradient correction

    real_t *K_traj;
    real_t *DK_traj;
    real_t *mu_traj;

    real_t *x;
    real_t *u;

    real_t *delta_DK_traj;
    real_t *adj_traj;
    real_t **jac_traj;

    real_t **sys_mat2;
    int_t **ipiv2;
    real_t **sys_sol2;

    struct blasfeo_dmat **str_mat2;
    struct blasfeo_dmat **str_sol2;

    // Memory for functions
    void *forward_vde_mem;
    void *jacobian_ode_mem;
} sim_lifted_irk_integrator_memory;



typedef struct {
    real_t *rhs_in;
    real_t *jac_tmp;
    real_t **VDE_tmp;
    real_t *out_tmp;
    int_t *ipiv;

    real_t *sys_mat;
    real_t *sys_sol;
    real_t *sys_sol_trans;

    real_t *trans;
    struct blasfeo_dmat *str_mat;
    struct blasfeo_dmat *str_sol;

    real_t *out_adj_tmp;

    // Workspace for functions
    void *forward_vde_work;
    void *jacobian_ode_work;
} sim_lifted_irk_integrator_workspace;



//
int sim_lifted_irk_integrator_calculate_args_size(sim_dims *dims, void *submodules_);
//
void *sim_lifted_irk_integrator_assign_args(sim_dims *dims, void **submodules_, void *raw_memory);
//
void *sim_lifted_irk_integrator_copy_args(sim_dims *dims, void *raw_memory, void *source_);
//
void sim_lifted_irk_integrator_initialize_default_args(sim_dims *dims, void *args_);
//
int sim_lifted_irk_integrator_calculate_memory_size(sim_dims *dims, void *args_);
//
void *sim_lifted_irk_integrator_assign_memory(sim_dims *dims, void *args_, void *raw_memory);
//
int sim_lifted_irk_integrator_calculate_workspace_size(sim_dims *dims, void *args_);
//
int sim_lifted_irk_integrator(sim_in *in, sim_out *out, void *args_, void *mem_, void *work_);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_LIFTED_IRK_INTEGRATOR_H_
