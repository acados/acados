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
#include "acados/sim/sim_common.h"
#include "acados/utils/external_function_generic.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/sim_interface.h"

// example specific
#include "{{ model.name }}_model/{{ model.name }}_model.h"
#include "acados_sim_solver_{{ model.name }}.h"

#define NX_    {{ dims.nx }}
#define NZ_    {{ dims.nz }}
#define NU_    {{ dims.nu }}
#define NP_    {{ dims.np }}

#define NX   NX_
#define NZ   NZ_
#define NU   NU_
#define NP   NP_

sim_config  * {{model.name}}_sim_config;
sim_in      * {{model.name}}_sim_in;
sim_out     * {{model.name}}_sim_out;
void        * {{model.name}}_sim_dims;
sim_opts    * {{model.name}}_sim_opts;
sim_solver  * {{model.name}}_sim_solver;

{% if solver_options.integrator_type == "ERK" %}
external_function_param_casadi * sim_forw_vde_casadi;
external_function_param_casadi * sim_expl_ode_fun_casadi;
{% elif solver_options.integrator_type == "IRK" %}
external_function_param_casadi * sim_impl_dae_fun;
external_function_param_casadi * sim_impl_dae_fun_jac_x_xdot_z;
external_function_param_casadi * sim_impl_dae_jac_x_xdot_u_z;
{% endif %}


int {{ model.name }}_acados_sim_create() {

    // initialize
    int ii;
    int jj;

    int nx = NX;
    int nu = NU;
    int nz = NZ;

    double Td = {{ solver_options.tf }}/ {{ dims.N }};

    {% if solver_options.integrator_type == "IRK" %}
    sim_impl_dae_fun = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi));
    sim_impl_dae_fun_jac_x_xdot_z = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi));
    sim_impl_dae_jac_x_xdot_u_z = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi));

    // external functions (implicit model)
    sim_impl_dae_fun->casadi_fun  = &{{ model.name }}_impl_dae_fun;
    sim_impl_dae_fun->casadi_work = &{{ model.name }}_impl_dae_fun_work;
    sim_impl_dae_fun->casadi_sparsity_in = &{{ model.name }}_impl_dae_fun_sparsity_in;
    sim_impl_dae_fun->casadi_sparsity_out = &{{ model.name }}_impl_dae_fun_sparsity_out;
    sim_impl_dae_fun->casadi_n_in = &{{ model.name }}_impl_dae_fun_n_in;
    sim_impl_dae_fun->casadi_n_out = &{{ model.name }}_impl_dae_fun_n_out;

    external_function_param_casadi_create(sim_impl_dae_fun, {{ dims.np }});

    sim_impl_dae_fun_jac_x_xdot_z->casadi_fun = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z;
    sim_impl_dae_fun_jac_x_xdot_z->casadi_work = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z_work;
    sim_impl_dae_fun_jac_x_xdot_z->casadi_sparsity_in = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z_sparsity_in;
    sim_impl_dae_fun_jac_x_xdot_z->casadi_sparsity_out = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z_sparsity_out;
    sim_impl_dae_fun_jac_x_xdot_z->casadi_n_in = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z_n_in;
    sim_impl_dae_fun_jac_x_xdot_z->casadi_n_out = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z_n_out;

    external_function_param_casadi_create(sim_impl_dae_fun_jac_x_xdot_z, {{ dims.np }});

    // external_function_param_casadi impl_dae_jac_x_xdot_u_z;
    sim_impl_dae_jac_x_xdot_u_z->casadi_fun = &{{ model.name }}_impl_dae_jac_x_xdot_u_z;
    sim_impl_dae_jac_x_xdot_u_z->casadi_work = &{{ model.name }}_impl_dae_jac_x_xdot_u_z_work;
    sim_impl_dae_jac_x_xdot_u_z->casadi_sparsity_in = &{{ model.name }}_impl_dae_jac_x_xdot_u_z_sparsity_in;
    sim_impl_dae_jac_x_xdot_u_z->casadi_sparsity_out = &{{ model.name }}_impl_dae_jac_x_xdot_u_z_sparsity_out;
    sim_impl_dae_jac_x_xdot_u_z->casadi_n_in = &{{ model.name }}_impl_dae_jac_x_xdot_u_z_n_in;
    sim_impl_dae_jac_x_xdot_u_z->casadi_n_out = &{{ model.name }}_impl_dae_jac_x_xdot_u_z_n_out;

    external_function_param_casadi_create(sim_impl_dae_jac_x_xdot_u_z, {{ dims.np }});

    {% elif solver_options.integrator_type == "ERK" %}
    // explicit ode
    sim_forw_vde_casadi = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi));
    sim_expl_ode_fun_casadi = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi));

    sim_forw_vde_casadi->casadi_fun = &{{ model.name }}_expl_vde_forw;
    sim_forw_vde_casadi->casadi_n_in = &{{ model.name }}_expl_vde_forw_n_in;
    sim_forw_vde_casadi->casadi_n_out = &{{ model.name }}_expl_vde_forw_n_out;
    sim_forw_vde_casadi->casadi_sparsity_in = &{{ model.name }}_expl_vde_forw_sparsity_in;
    sim_forw_vde_casadi->casadi_sparsity_out = &{{ model.name }}_expl_vde_forw_sparsity_out;
    sim_forw_vde_casadi->casadi_work = &{{ model.name }}_expl_vde_forw_work;

    external_function_param_casadi_create(sim_forw_vde_casadi, {{ dims.np }});

    sim_expl_ode_fun_casadi->casadi_fun = &{{ model.name }}_expl_ode_fun;
    sim_expl_ode_fun_casadi->casadi_n_in = &{{ model.name }}_expl_ode_fun_n_in;
    sim_expl_ode_fun_casadi->casadi_n_out = &{{ model.name }}_expl_ode_fun_n_out;
    sim_expl_ode_fun_casadi->casadi_sparsity_in = &{{ model.name }}_expl_ode_fun_sparsity_in;
    sim_expl_ode_fun_casadi->casadi_sparsity_out = &{{ model.name }}_expl_ode_fun_sparsity_out;
    sim_expl_ode_fun_casadi->casadi_work = &{{ model.name }}_expl_ode_fun_work;

    external_function_param_casadi_create(sim_expl_ode_fun_casadi, {{ dims.np }});

    {% endif %}

    // sim plan & config

    // choose plan
    sim_solver_plan plan;

    plan.sim_solver = {{ solver_options.integrator_type }};

    // create correct config based on plan
    {{ model.name }}_sim_config = sim_config_create(plan);

    // sim dims

    {{ model.name }}_sim_dims = sim_dims_create({{ model.name }}_sim_config);
    sim_dims_set({{ model.name }}_sim_config, {{ model.name }}_sim_dims, "nx", &nx);
    sim_dims_set({{ model.name }}_sim_config, {{ model.name }}_sim_dims, "nu", &nu);
    sim_dims_set({{ model.name }}_sim_config, {{ model.name }}_sim_dims, "nz", &nz);

    // sim opts

    {{ model.name }}_sim_opts = sim_opts_create({{ model.name }}_sim_config, {{ model.name }}_sim_dims);

    {{ model.name }}_sim_opts->ns = {{ solver_options.sim_method_num_stages }}; // number of stages in rk integrator
    {{ model.name }}_sim_opts->num_steps = {{ solver_options.sim_method_num_steps }}; // number of integration steps
    {{ model.name }}_sim_opts->sens_adj = false;
    {{ model.name }}_sim_opts->sens_forw = true;
    {% if solver_options.integrator_type == "IRK" %}
    {{ model.name }}_sim_opts->sens_algebraic = false;
    {{ model.name }}_sim_opts->output_z = false;
    {{ model.name }}_sim_opts->newton_iter = 5;
    {{ model.name }}_sim_opts->jac_reuse = false;
    {% endif %}

    // sim in / out

    {{ model.name }}_sim_in  = sim_in_create({{ model.name }}_sim_config, {{ model.name }}_sim_dims);
    {{ model.name }}_sim_out = sim_out_create({{ model.name }}_sim_config, {{ model.name }}_sim_dims);

    {{ model.name }}_sim_in->T = Td;

    // external functions
    {% if solver_options.integrator_type == "IRK" %}
    {{ model.name }}_sim_config->model_set({{ model.name }}_sim_in->model, "impl_ode_fun", sim_impl_dae_fun);
    {{ model.name }}_sim_config->model_set({{ model.name }}_sim_in->model, "impl_ode_fun_jac_x_xdot", sim_impl_dae_fun_jac_x_xdot_z);
    {{ model.name }}_sim_config->model_set({{ model.name }}_sim_in->model, "impl_ode_jac_x_xdot_u", sim_impl_dae_jac_x_xdot_u_z);
    {% elif solver_options.integrator_type == "ERK" %}
    {{ model.name }}_sim_config->model_set({{ model.name }}_sim_in->model, "expl_vde_for", sim_forw_vde_casadi);
    {{ model.name }}_sim_config->model_set({{ model.name }}_sim_in->model, "expl_ode_fun", sim_expl_ode_fun_casadi);
    {% endif %}

    // sim solver

    {{ model.name }}_sim_solver = sim_solver_create({{ model.name }}_sim_config, {{ model.name }}_sim_dims, {{ model.name }}_sim_opts);

    // initialize state and input to zero
    // x
    for (ii = 0; ii < NX; ii++)
        {{ model.name }}_sim_in->x[ii] = 0.0;

    // u
    for (ii = 0; ii < NU; ii++)
        {{ model.name }}_sim_in->u[ii] = 0.0;
    int status = 0;

    return status;
}


int {{ model.name }}_acados_sim_solve() {

    // integrate dynamics using acados sim_solver
    int status = sim_solve({{ model.name }}_sim_solver, {{ model.name }}_sim_in, {{ model.name }}_sim_out);
    if (status != 0)
        printf("error in {{ model.name }}_acados_sim_solve()! Exiting.\n");

    return status;
}


int {{ model.name }}_acados_sim_free() {

    // free memory
    sim_solver_destroy({{ model.name }}_sim_solver);
    sim_in_destroy({{ model.name }}_sim_in);
    sim_out_destroy({{ model.name }}_sim_out);
    sim_opts_destroy({{ model.name }}_sim_opts);
    sim_dims_destroy({{ model.name }}_sim_dims);
    sim_config_destroy({{ model.name }}_sim_config);

    // free external function
    {% if solver_options.integrator_type == "IRK" %}
        external_function_param_casadi_free(sim_impl_dae_fun);
        external_function_param_casadi_free(sim_impl_dae_fun_jac_x_xdot_z);
        external_function_param_casadi_free(sim_impl_dae_jac_x_xdot_u_z);
    {% elif solver_options.integrator_type == "ERK" %}
        external_function_param_casadi_free(sim_forw_vde_casadi);
        external_function_param_casadi_free(sim_expl_ode_fun_casadi);
    {% endif %}

    return 0;
}


/* getters pointers to C objects*/
sim_config * {{ model.name }}_acados_get_sim_config()
{
    return {{ model.name }}_sim_config;
};

sim_in * {{ model.name }}_acados_get_sim_in()
{
    return {{ model.name }}_sim_in;
};

sim_out * {{ model.name }}_acados_get_sim_out()
{
    return {{ model.name }}_sim_out;
};

void * {{ model.name }}_acados_get_sim_dims()
{
    return {{ model.name }}_sim_dims;
};

sim_opts * {{ model.name }}_acados_get_sim_opts()
{
    return {{ model.name }}_sim_opts;
};

sim_solver  * {{ model.name }}_acados_get_sim_solver()
{
    return {{ model.name }}_sim_solver;
};

