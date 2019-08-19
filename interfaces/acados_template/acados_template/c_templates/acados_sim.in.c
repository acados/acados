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

// standard
#include <stdio.h>
#include <stdlib.h>
// acados
#include "acados/utils/print.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// example specific

#include "acados_sim_{{ ocp.model_name }}.h"

#define NX_    {{ ocp.dims.nx }}
#define NZ_    {{ ocp.dims.nz }}
#define NU_    {{ ocp.dims.nu }}
#define NP_    {{ ocp.dims.np }}

#if NX_ < 1
#define NX   1
#else
#define NX   NX_
#endif

#if NZ_ < 1
#define NZ   1
#else
#define NZ   NZ_
#endif

#if NU_ < 1
#define NU   1
#else
#define NU   NU_
#endif

#if NP_ < 1
#define NP   1
#else
#define NP   NP_
#endif

int {{ ocp.model_name}}_acados_sim_create() {

	// initialize

    int ii;
    int jj;

    int nx = NX;
    int nu = NU;
    int nz = NZ;

    double Td = {{ ocp.solver_config.tf }}/N;

	// external functions (implicit model)
	sim_impl_dae_fun.casadi_fun  = &{{ ocp.model_name }}_impl_dae_fun;
	sim_impl_dae_fun.casadi_work = &{{ ocp.model_name }}_impl_dae_fun_work;
	sim_impl_dae_fun.casadi_sparsity_in = &{{ ocp.model_name }}_impl_dae_fun_sparsity_in;
	sim_impl_dae_fun.casadi_sparsity_out = &{{ ocp.model_name }}_impl_dae_fun_sparsity_out;
	sim_impl_dae_fun.casadi_n_in = &{{ ocp.model_name }}_impl_dae_fun_n_in;
	sim_impl_dae_fun.casadi_n_out = &{{ ocp.model_name }}_impl_dae_fun_n_out;
	external_function_param_casadi_create(&sim_impl_dae_fun, 3);

	// external_function_casadi impl_dae_fun_jac_x_xdot_z;
	sim_impl_dae_fun_jac_x_xdot_z.casadi_fun = &{{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_z;
	sim_impl_dae_fun_jac_x_xdot_z.casadi_work = &{{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_z_work;
	sim_impl_dae_fun_jac_x_xdot_z.casadi_sparsity_in = &{{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_z_sparsity_in;
	sim_impl_dae_fun_jac_x_xdot_z.casadi_sparsity_out = &{{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_z_sparsity_out;
	sim_impl_dae_fun_jac_x_xdot_z.casadi_n_in = &{{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_z_n_in;
	sim_impl_dae_fun_jac_x_xdot_z.casadi_n_out = &{{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_z_n_out;
	external_function_param_casadi_create(&sim_impl_dae_fun_jac_x_xdot_z, 3);

	// external_function_casadi impl_dae_jac_x_xdot_u_z;
	sim_impl_dae_jac_x_xdot_u_z.casadi_fun = &{{ ocp.model_name }}_impl_dae_jac_x_xdot_u_z;
	sim_impl_dae_jac_x_xdot_u_z.casadi_work = &{{ ocp.model_name }}_impl_dae_jac_x_xdot_u_z_work;
	sim_impl_dae_jac_x_xdot_u_z.casadi_sparsity_in = &{{ ocp.model_name }}_impl_dae_jac_x_xdot_u_z_sparsity_in;
	sim_impl_dae_jac_x_xdot_u_z.casadi_sparsity_out = &{{ ocp.model_name }}_impl_dae_jac_x_xdot_u_z_sparsity_out;
	sim_impl_dae_jac_x_xdot_u_z.casadi_n_in = &{{ ocp.model_name }}_impl_dae_jac_x_xdot_u_z_n_in;
	sim_impl_dae_jac_x_xdot_u_z.casadi_n_out = &{{ ocp.model_name }}_impl_dae_jac_x_xdot_u_z_n_out;
	external_function_param_casadi_create(&sim_impl_dae_jac_x_xdot_u_z, 3);

    // sim plan & config

    // choose plan
    sim_solver_plan plan;

    plan.sim_solver = {{ ocp.solver_config.integrator_type }};

    // create correct config based on plan
    {{ ocp.model_name }}_sim_config = sim_config_create(plan);

    // sim dims

    {{ ocp.model_name }}_sim_dims = sim_dims_create({{ ocp.model_name }}_sim_config);
    sim_dims_set({{ ocp.model_name }}_sim_config, {{ ocp.model_name }}_sim_dims, "nx", &nx);
    sim_dims_set({{ ocp.model_name }}_sim_config, {{ ocp.model_name }}_sim_dims, "nu", &nu);
    sim_dims_set({{ ocp.model_name }}_sim_config, {{ ocp.model_name }}_sim_dims, "nz", &nz);

    // sim opts

    {{ ocp.model_name }}_sim_opts = sim_opts_create({{ ocp.model_name }}_sim_config, {{ ocp.model_name }}_sim_dims);

    {{ ocp.model_name }}_sim_opts->ns = 1; // number of stages in rk integrator
    {{ ocp.model_name }}_sim_opts->num_steps = 1; // number of integration steps
    {% if ocp.solver_config.integrator_type == "IRK" %}
    {{ ocp.model_name }}_sim_opts->sens_algebraic = false;
    {{ ocp.model_name }}_sim_opts->sens_adj = false;
    {{ ocp.model_name }}_sim_opts->sens_forw = false;
    {{ ocp.model_name }}_sim_opts->output_z = false;
    {{ ocp.model_name }}_sim_opts->newton_iter = 5;
    {{ ocp.model_name }}_sim_opts->jac_reuse = false;
    {% endif %}

    // sim in / out

    {{ ocp.model_name }}_sim_in  = sim_in_create({{ ocp.model_name }}_sim_config, {{ ocp.model_name }}_sim_dims);
    {{ ocp.model_name }}_sim_out = sim_out_create({{ ocp.model_name }}_sim_config, {{ ocp.model_name }}_sim_dims);

    {{ ocp.model_name }}_sim_in->T = Td;

    // external functions
    sim_model_set({{ ocp.model_name }}_sim_config, {{ ocp.model_name }}_sim_in, "impl_ode_fun", &sim_impl_dae_fun);
    sim_model_set({{ ocp.model_name }}_sim_config, {{ ocp.model_name }}_sim_in, "impl_ode_fun_jac_x_xdot", &sim_impl_dae_fun_jac_x_xdot_z);
    sim_model_set({{ ocp.model_name }}_sim_config, {{ ocp.model_name }}_sim_in, "impl_ode_jac_x_xdot_u", &sim_impl_dae_jac_x_xdot_u_z);

    // sim solver

    {{ ocp.model_name }}_sim_solver = sim_create({{ ocp.model_name }}_sim_config, {{ ocp.model_name }}_sim_dims, {{ ocp.model_name }}_sim_opts);
    
    // initialize state and input to zero
    // x
    for (ii = 0; ii < NX; ii++)
        {{ ocp.model_name }}_sim_in->x[ii] = 0.0;
    
    // u
    for (ii = 0; ii < NU; ii++)
        {{ ocp.model_name }}_sim_in->u[ii] = 0.0;
    int status = 0;

    return status;
}

int {{ ocp.model_name }}_acados_sim_solve() {

    // integrate dynamics using acados sim_solver 
    int status = sim_solve({{ ocp.model_name }}_sim_solver, {{ ocp.model_name }}_sim_in, {{ ocp.model_name }}_sim_out);
    if (status != 0)
        printf("error in {{ ocp.model_name }}_acados_sim_solve()! Exiting.\n");

    return status;
}

int {{ ocp.model_name }}_acados_sim_free() {

    // free memory
    sim_solver_destroy({{ ocp.model_name }}_sim_solver);
    sim_in_destroy({{ ocp.model_name }}_sim_in);
    sim_out_destroy({{ ocp.model_name }}_sim_out);
    sim_opts_destroy({{ ocp.model_name }}_sim_opts);
    sim_dims_destroy({{ ocp.model_name }}_sim_dims);
    sim_config_destroy({{ ocp.model_name }}_sim_config);

    // free external function 
    {% if ocp.solver_config.integrator_type == "IRK" %}
    for(int i = 0; i < {{ocp.dims.N}}; i++) {
        {% if ocp.dims.np < 1 %}
        external_function_casadi_free(&sim_impl_dae_fun[i]);
        external_function_casadi_free(&sim_impl_dae_fun_jac_x_xdot_z[i]);
        external_function_casadi_free(&sim_impl_dae_jac_x_xdot_u_z[i]);
        {% else %}
        external_function_param_casadi_free(&sim_impl_dae_fun[i]);
        external_function_param_casadi_free(&sim_impl_dae_fun_jac_x_xdot_z[i]);
        external_function_param_casadi_free(&sim_impl_dae_jac_x_xdot_u_z[i]);
        {% endif %}
    }
    {% else %}
    for(int i = 0; i < {{ocp.dims.N}}; i++) {
        {% if ocp.dims.np < 1 %}
        external_function_casadi_free(&sim_forw_vde_casadi[i]);
        {% else %}
        external_function_param_casadi_free(&sim_forw_vde_casadi[i]);
        {% endif %}
    }
    {% endif %}
    
    return 0;
}

sim_config  * {{ocp.model_name}}_acados_get_sim_config() {  return {{ocp.model_name}}_sim__config()};
sim_in      * {{ocp.model_name}}_acados_get_sim_in(){       return {{ocp.model_name}}_sim__in()};
sim_out     * {{ocp.model_name}}_acados_get_sim_out(){      return {{ocp.model_name}}_sim__out()};
void        * {{ocp.model_name}}_acados_get_sim_dims(){     return {{ocp.model_name}}_sim__dims()};
sim_opts    * {{ocp.model_name}}_acados_get_sim_opts{       return {{ocp.model_name}}_sim__opts};
sim_solver  * {{ocp.model_name}}_acados_get_sim_solver(){   return {{ocp.model_name}}_sim__solver()};

