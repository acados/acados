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
#include "acados/utils/types.h"

typedef struct {

} sim_irk_memory;

typedef struct {

    struct d_strmat *JG; // jacobian of G over K (nx*ns, nx*ns)
    struct d_strmat *JGf; // jacobian of G over x and u (nx*ns, nx+nu);
    struct d_strmat *JKf; // jacobian of K over x and u (nx*ns, nx+nu);
    struct d_strmat *JFK; // jacobian of F over K (nx, nx*ns) 
    struct d_strmat *S_forw; // forward sensitivities

    struct d_strvec *rG; // residuals of G (nx*ns)
    struct d_strvec *K; // internal variables (nx*ns)
    struct d_strvec *xt; // temporary x
    struct d_strvec *xn; // x at each integration step
    
    struct d_strvec *lambda; // adjoint seed (nx+nu)
    struct d_strvec *lambdaK; // auxiliary variable (nx*ns)
    
    double *rGt; // temporary residuals of G (nx, 1)
    double *jac_out; // temporary Jacobian of ode (nx, 2*nx+nu)
    double *Jt; // temporary Jacobian of ode (nx, nx)
    double *ode_args; // pointer to ode args
    double *S_adj_w;
    int *ipiv; // index of pivot vector

    struct d_strvec **xn_traj; // xn trajectory
    struct d_strvec **K_traj;  // K trajectory
    struct d_strmat **JG_traj; // JG trajectory

} sim_irk_workspace;

int sim_irk_opts_calculate_size(sim_dims *dims);

void *sim_irk_assign_opts(sim_dims *dims, void *raw_memory);

void sim_irk_initialize_default_args(sim_dims *dims, void *opts_);

int sim_irk_calculate_memory_size(sim_dims *dims, void *opts_);

void *sim_irk_assign_memory(sim_dims *dims, void *opts_, void *raw_memory);

void *sim_irk_create_memory(sim_dims *dims, void *opts_);

int sim_irk(sim_in *in, sim_out *out, void *opts_, void *mem_, void *work_);

int sim_irk_calculate_workspace_size(sim_dims *dims, void *opts_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_IRK_INTEGRATOR_H_
