/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */

// standard
#include <stdio.h>
#include <stdlib.h>
// acados
#include "acados/utils/print.h"
#include "acados_c/sim_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/external_function_generic.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/sim_interface.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// example specific

#include "{{ ocp.model.name }}_model/{{ ocp.model.name }}_model.h"
#include "acados_sim_solver_{{ ocp.model.name }}.h"

#define NX_    {{ ocp.dims.nx }}
#define NZ_    {{ ocp.dims.nz }}
#define NU_    {{ ocp.dims.nu }}
#define NP_    {{ ocp.dims.np }}

#define NX   NX_
#define NZ   NZ_
#define NU   NU_
#define NP   NP_

int {{ ocp.model.name}}_acados_sim_create() {

	// initialize

    int ii;
    int jj;

    int nx = NX;
    int nu = NU;
    int nz = NZ;

    double Td = {{ ocp.solver_options.tf }}/ {{ ocp.dims.N }};

    {% if ocp.solver_options.integrator_type == "IRK" %}

    sim_impl_dae_fun = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi));
    sim_impl_dae_fun_jac_x_xdot_z = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi));
    sim_impl_dae_jac_x_xdot_u_z = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi));

	// external functions (implicit model)
	sim_impl_dae_fun->casadi_fun  = &{{ ocp.model.name }}_impl_dae_fun;
	sim_impl_dae_fun->casadi_work = &{{ ocp.model.name }}_impl_dae_fun_work;
	sim_impl_dae_fun->casadi_sparsity_in = &{{ ocp.model.name }}_impl_dae_fun_sparsity_in;
	sim_impl_dae_fun->casadi_sparsity_out = &{{ ocp.model.name }}_impl_dae_fun_sparsity_out;
	sim_impl_dae_fun->casadi_n_in = &{{ ocp.model.name }}_impl_dae_fun_n_in;
	sim_impl_dae_fun->casadi_n_out = &{{ ocp.model.name }}_impl_dae_fun_n_out;

    external_function_param_casadi_create(sim_impl_dae_fun, {{ ocp.dims.np }});

	// external_function_param_casadi impl_dae_fun_jac_x_xdot_z;
	sim_impl_dae_fun_jac_x_xdot_z->casadi_fun = &{{ ocp.model.name }}_impl_dae_fun_jac_x_xdot_z;
	sim_impl_dae_fun_jac_x_xdot_z->casadi_work = &{{ ocp.model.name }}_impl_dae_fun_jac_x_xdot_z_work;
	sim_impl_dae_fun_jac_x_xdot_z->casadi_sparsity_in = &{{ ocp.model.name }}_impl_dae_fun_jac_x_xdot_z_sparsity_in;
	sim_impl_dae_fun_jac_x_xdot_z->casadi_sparsity_out = &{{ ocp.model.name }}_impl_dae_fun_jac_x_xdot_z_sparsity_out;
	sim_impl_dae_fun_jac_x_xdot_z->casadi_n_in = &{{ ocp.model.name }}_impl_dae_fun_jac_x_xdot_z_n_in;
	sim_impl_dae_fun_jac_x_xdot_z->casadi_n_out = &{{ ocp.model.name }}_impl_dae_fun_jac_x_xdot_z_n_out;

    external_function_param_casadi_create(sim_impl_dae_fun_jac_x_xdot_z, {{ ocp.dims.np }});

	// external_function_param_casadi impl_dae_jac_x_xdot_u_z;
	sim_impl_dae_jac_x_xdot_u_z->casadi_fun = &{{ ocp.model.name }}_impl_dae_jac_x_xdot_u_z;
	sim_impl_dae_jac_x_xdot_u_z->casadi_work = &{{ ocp.model.name }}_impl_dae_jac_x_xdot_u_z_work;
	sim_impl_dae_jac_x_xdot_u_z->casadi_sparsity_in = &{{ ocp.model.name }}_impl_dae_jac_x_xdot_u_z_sparsity_in;
	sim_impl_dae_jac_x_xdot_u_z->casadi_sparsity_out = &{{ ocp.model.name }}_impl_dae_jac_x_xdot_u_z_sparsity_out;
	sim_impl_dae_jac_x_xdot_u_z->casadi_n_in = &{{ ocp.model.name }}_impl_dae_jac_x_xdot_u_z_n_in;
	sim_impl_dae_jac_x_xdot_u_z->casadi_n_out = &{{ ocp.model.name }}_impl_dae_jac_x_xdot_u_z_n_out;

    external_function_param_casadi_create(sim_impl_dae_jac_x_xdot_u_z, {{ ocp.dims.np }});

    {% else %}
    // explicit ode
    sim_forw_vde_casadi = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi));
    sim_expl_ode_fun_casadi = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi));

    sim_forw_vde_casadi->casadi_fun = &{{ ocp.model.name }}_expl_vde_forw;
    sim_forw_vde_casadi->casadi_n_in = &{{ ocp.model.name }}_expl_vde_forw_n_in;
    sim_forw_vde_casadi->casadi_n_out = &{{ ocp.model.name }}_expl_vde_forw_n_out;
    sim_forw_vde_casadi->casadi_sparsity_in = &{{ ocp.model.name }}_expl_vde_forw_sparsity_in;
    sim_forw_vde_casadi->casadi_sparsity_out = &{{ ocp.model.name }}_expl_vde_forw_sparsity_out;
    sim_forw_vde_casadi->casadi_work = &{{ ocp.model.name }}_expl_vde_forw_work;

    external_function_param_casadi_create(sim_forw_vde_casadi, {{ ocp.dims.np }});

    sim_expl_ode_fun_casadi->casadi_fun = &{{ ocp.model.name }}_expl_ode_fun;
    sim_expl_ode_fun_casadi->casadi_n_in = &{{ ocp.model.name }}_expl_ode_fun_n_in;
    sim_expl_ode_fun_casadi->casadi_n_out = &{{ ocp.model.name }}_expl_ode_fun_n_out;
    sim_expl_ode_fun_casadi->casadi_sparsity_in = &{{ ocp.model.name }}_expl_ode_fun_sparsity_in;
    sim_expl_ode_fun_casadi->casadi_sparsity_out = &{{ ocp.model.name }}_expl_ode_fun_sparsity_out;
    sim_expl_ode_fun_casadi->casadi_work = &{{ ocp.model.name }}_expl_ode_fun_work;

    external_function_param_casadi_create(sim_expl_ode_fun_casadi, {{ ocp.dims.np }});

    {% endif %}

    // sim plan & config

    // choose plan
    sim_solver_plan_t plan;

    plan.sim_solver = {{ ocp.solver_options.integrator_type }};

    // create correct config based on plan
    {{ ocp.model.name }}_sim_config = sim_config_create(plan);

    // sim dims

    {{ ocp.model.name }}_sim_dims = sim_dims_create({{ ocp.model.name }}_sim_config);
    sim_dims_set({{ ocp.model.name }}_sim_config, {{ ocp.model.name }}_sim_dims, "nx", &nx);
    sim_dims_set({{ ocp.model.name }}_sim_config, {{ ocp.model.name }}_sim_dims, "nu", &nu);
    sim_dims_set({{ ocp.model.name }}_sim_config, {{ ocp.model.name }}_sim_dims, "nz", &nz);

    // sim opts

    {{ ocp.model.name }}_sim_opts = sim_opts_create({{ ocp.model.name }}_sim_config, {{ ocp.model.name }}_sim_dims);

    {{ ocp.model.name }}_sim_opts->ns = {{ ocp.solver_options.sim_method_num_stages }}; // number of stages in rk integrator
    {{ ocp.model.name }}_sim_opts->num_steps = {{ ocp.solver_options.sim_method_num_steps }}; // number of integration steps
    {{ ocp.model.name }}_sim_opts->sens_adj = false;
    {{ ocp.model.name }}_sim_opts->sens_forw = true;
    {% if ocp.solver_options.integrator_type == "IRK" %}
    {{ ocp.model.name }}_sim_opts->sens_algebraic = false;
    {{ ocp.model.name }}_sim_opts->output_z = false;
    {{ ocp.model.name }}_sim_opts->newton_iter = 5;
    {{ ocp.model.name }}_sim_opts->jac_reuse = false;
    {% endif %}

    // sim in / out

    {{ ocp.model.name }}_sim_in  = sim_in_create({{ ocp.model.name }}_sim_config, {{ ocp.model.name }}_sim_dims);
    {{ ocp.model.name }}_sim_out = sim_out_create({{ ocp.model.name }}_sim_config, {{ ocp.model.name }}_sim_dims);

    {{ ocp.model.name }}_sim_in->T = Td;

    // external functions
    {% if ocp.solver_options.integrator_type == "IRK" %}
    {{ ocp.model.name }}_sim_config->model_set({{ ocp.model.name }}_sim_in->model, "impl_ode_fun", sim_impl_dae_fun);
    {{ ocp.model.name }}_sim_config->model_set({{ ocp.model.name }}_sim_in->model, "impl_ode_fun_jac_x_xdot", sim_impl_dae_fun_jac_x_xdot_z);
    {{ ocp.model.name }}_sim_config->model_set({{ ocp.model.name }}_sim_in->model, "impl_ode_jac_x_xdot_u", sim_impl_dae_jac_x_xdot_u_z);

    {% else %}
    {% if ocp.solver_options.integrator_type == "ERK" %} 
    {{ ocp.model.name }}_sim_config->model_set({{ ocp.model.name }}_sim_in->model, "expl_vde_for", sim_forw_vde_casadi);
    {{ ocp.model.name }}_sim_config->model_set({{ ocp.model.name }}_sim_in->model, "expl_ode_fun", sim_expl_ode_fun_casadi);
    {% endif %}
    {% endif %}

    // sim solver

    {{ ocp.model.name }}_sim_solver = sim_solver_create({{ ocp.model.name }}_sim_config, {{ ocp.model.name }}_sim_dims, {{ ocp.model.name }}_sim_opts);
    
    // initialize state and input to zero
    // x
    for (ii = 0; ii < NX; ii++)
        {{ ocp.model.name }}_sim_in->x[ii] = 0.0;
    
    // u
    for (ii = 0; ii < NU; ii++)
        {{ ocp.model.name }}_sim_in->u[ii] = 0.0;
    int status = 0;

    return status;
}

int {{ ocp.model.name }}_acados_sim_solve() {

    // integrate dynamics using acados sim_solver 
    int status = sim_solve({{ ocp.model.name }}_sim_solver, {{ ocp.model.name }}_sim_in, {{ ocp.model.name }}_sim_out);
    if (status != 0)
        printf("error in {{ ocp.model.name }}_acados_sim_solve()! Exiting.\n");

    return status;
}

int {{ ocp.model.name }}_acados_sim_free() {

    // free memory
    sim_solver_destroy({{ ocp.model.name }}_sim_solver);
    sim_in_destroy({{ ocp.model.name }}_sim_in);
    sim_out_destroy({{ ocp.model.name }}_sim_out);
    sim_opts_destroy({{ ocp.model.name }}_sim_opts);
    sim_dims_destroy({{ ocp.model.name }}_sim_dims);
    sim_config_destroy({{ ocp.model.name }}_sim_config);

    // free external function 
    {% if ocp.solver_options.integrator_type == "IRK" %}
        external_function_param_casadi_free(sim_impl_dae_fun);
        external_function_param_casadi_free(sim_impl_dae_fun_jac_x_xdot_z);
        external_function_param_casadi_free(sim_impl_dae_jac_x_xdot_u_z);
    {% else %}
        external_function_param_casadi_free(sim_forw_vde_casadi);
        external_function_param_casadi_free(sim_expl_ode_fun_casadi);
    {% endif %}
    
    return 0;
}

sim_config  * {{ocp.model.name}}_acados_get_sim_config() {  
    return {{ocp.model.name}}_sim_config; };

sim_in      * {{ocp.model.name}}_acados_get_sim_in(){       
    return {{ocp.model.name}}_sim_in; };

sim_out     * {{ocp.model.name}}_acados_get_sim_out(){      
    return {{ocp.model.name}}_sim_out; };

void        * {{ocp.model.name}}_acados_get_sim_dims(){     
    return {{ocp.model.name}}_sim_dims; };

sim_opts    * {{ocp.model.name}}_acados_get_sim_opts(){     
    return {{ocp.model.name}}_sim_opts; };

sim_solver  * {{ocp.model.name}}_acados_get_sim_solver(){   
    return {{ocp.model.name}}_sim_solver; };

