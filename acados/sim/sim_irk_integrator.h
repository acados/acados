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

#ifndef ACADOS_SIM_SIM_IRK_INTEGRATOR_H_
#define ACADOS_SIM_SIM_IRK_INTEGRATOR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/sim/sim_common.h"
#include "acados/utils/external_function.h"
#include "acados/utils/types.h"



typedef struct {
    external_function_fcn_ptrs *impl_res;
    external_function_fcn_ptrs *impl_jac;
} sim_irk_integrator_submodules;



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
    bool jac_reuse;

    // Submodules
    sim_irk_integrator_submodules submodules;

    // Arguments for functions
    void *impl_res_args;
    void *impl_jac_args;
} sim_irk_integrator_args;



typedef struct {
    // Memory for functions
    void *impl_res_mem;
    void *impl_jac_mem;
} sim_irk_integrator_memory;



typedef struct {

    struct blasfeo_dmat *JGK; // jacobian of G over K (nx*ns, nx*ns)
    struct blasfeo_dmat *JGf; // jacobian of G over x and u (nx*ns, nx+nu);
    struct blasfeo_dmat *JKf; // jacobian of K over x and u (nx*ns, nx+nu);
    // struct blasfeo_dmat *JFK; // jacobian of F over K (nx, nx*ns) 
    struct blasfeo_dmat *S_forw; // forward sensitivities

    struct blasfeo_dvec *rG; // residuals of G (nx*ns)
    struct blasfeo_dvec *K; // internal variables (nx*ns)
    struct blasfeo_dvec *xt; // temporary x
    struct blasfeo_dvec *xn; // x at each integration step
    
    struct blasfeo_dvec *lambda; // adjoint seed (nx+nu)
    struct blasfeo_dvec *lambdaK; // auxiliary variable (nx*ns)
    
    double *rGt; // temporary residuals of G (nx, 1)
    double *jac_out; // temporary Jacobian of ode (nx, 2*nx+nu)
    double *Jt; // temporary Jacobian of ode (nx, nx)
    double *ode_args; // pointer to ode args
    double *S_adj_w;
    int *ipiv; // index of pivot vector

    struct blasfeo_dvec *xn_traj; // xn trajectory
    struct blasfeo_dvec *K_traj;  // K trajectory
    struct blasfeo_dmat *JG_traj; // JG trajectory

    // Workspace for functions
    void *impl_res_work;
    void *impl_jac_work;

} sim_irk_integrator_workspace;



//
int sim_irk_integrator_calculate_args_size(sim_dims *dims, void *submodules_);
//
void *sim_irk_integrator_assign_args(sim_dims *dims, void **submodules_, void *raw_memory);
//
void *sim_irk_integrator_copy_args(sim_dims *dims, void *raw_memory, void *source_);
//
void sim_irk_integrator_initialize_default_args(sim_dims *dims, void *args_);
//
int sim_irk_integrator_calculate_memory_size(sim_dims *dims, void *args_);
//
void *sim_irk_integrator_assign_memory(sim_dims *dims, void *args_, void *raw_memory);
//
int sim_irk_integrator_calculate_workspace_size(sim_dims *dims, void *args_);
//
int sim_irk_integrator(sim_in *in, sim_out *out, void *args_, void *mem_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_IRK_INTEGRATOR_H_
