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
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados_solver_{{ocp.model.name}}.h"
#include "acados_sim_solver_{{ocp.model.name}}.h"


// ** global data **
ocp_nlp_in * nlp_in;
ocp_nlp_out * nlp_out;
ocp_nlp_solver * nlp_solver;
void * nlp_opts;
ocp_nlp_plan * nlp_solver_plan;
ocp_nlp_config * nlp_config;
ocp_nlp_dims * nlp_dims;

sim_config  * {{ocp.model.name}}_sim_config;
sim_in      * {{ocp.model.name}}_sim_in;
sim_out     * {{ocp.model.name}}_sim_out; 
void        * {{ocp.model.name}}_sim_dims;
sim_opts    * {{ocp.model.name}}_sim_opts;
sim_solver  * {{ocp.model.name}}_sim_solver; 

{% if ocp.solver_config.integrator_type == "ERK" %}
{% if ocp.dims.np < 1 %}
external_function_param_casadi * forw_vde_casadi;
external_function_param_casadi * sim_forw_vde_casadi;
external_function_param_casadi * sim_expl_ode_fun_casadi;
{% else %}
external_function_param_casadi * forw_vde_casadi;
external_function_param_casadi * sim_forw_vde_casadi;
external_function_param_casadi * sim_expl_ode_fun_casadi;
{% endif %}
{% if ocp.solver_config.hessian_approx == "EXACT" %} 
{% if ocp.dims.np < 1 %}
external_function_param_casadi * hess_vde_casadi;
{% else %}
external_function_param_casadi * hess_vde_casadi;
{% endif %}
{% endif %}
{% else %}
{% if ocp.solver_config.integrator_type == "IRK" %}
{% if ocp.dims.np < 1 %}
external_function_param_casadi * impl_dae_fun;
external_function_param_casadi * impl_dae_fun_jac_x_xdot_z;
external_function_param_casadi * impl_dae_jac_x_xdot_u_z;
external_function_param_casadi * sim_impl_dae_fun;
external_function_param_casadi * sim_impl_dae_fun_jac_x_xdot_z;
external_function_param_casadi * sim_impl_dae_jac_x_xdot_u_z;
{% else %}
external_function_param_casadi * impl_dae_fun;
external_function_param_casadi * impl_dae_fun_jac_x_xdot_z;
external_function_param_casadi * impl_dae_jac_x_xdot_u_z;
external_function_param_casadi * sim_impl_dae_fun;
external_function_param_casadi * sim_impl_dae_fun_jac_x_xdot_z;
external_function_param_casadi * sim_impl_dae_jac_x_xdot_u_z;
{% endif %}
{% endif %}
{% endif %}
{% if ocp.dims.npd > 0 %}
external_function_param_casadi * p_constraint;
{% endif %}
{% if ocp.dims.npd_e > 0 %}
external_function_param_casadi p_e_constraint;
{% endif %}
{% if ocp.dims.nh > 0 %}
external_function_param_casadi * h_constraint;
{% endif %}
{% if ocp.dims.nh_e > 0 %}
external_function_param_casadi h_e_constraint;
{% endif %}

int main() {

    // test integrator first
    int sim_status = 0;
    sim_status = {{ ocp.model.name }}_acados_sim_create();

    // set sim input
    double x_sim[{{ocp.dims.nx}}];
    {% for item in ocp.constraints.x0 %}
    x_sim[{{ loop.index0 }}] = {{ item }};
    {% endfor %}

    // set initial condition
    double u_sim[{{ocp.dims.nu}}];
    {% for item in range(ocp.dims.nu) %}
    u_sim[{{ loop.index0 }}] = {{ item }};
    {% endfor %}
    
    // seeds forw
    for (int ii = 0; ii < {{ocp.dims.nx}} * ({{ocp.dims.nx}} + {{ocp.dims.nu}}); ii++)
        {{ ocp.model.name }}_sim_in->S_forw[ii] = 0.0;
    for (int ii = 0; ii < {{ocp.dims.nx}}; ii++)
        {{ ocp.model.name }}_sim_in->S_forw[ii * ({{ocp.dims.nx}} + 1)] = 1.0;

    double Td = {{ ocp.solver_config.tf }}/ {{ ocp.dims.N }};
    sim_in_set({{ ocp.model.name }}_sim_config, {{ ocp.model.name }}_sim_dims, {{ ocp.model.name }}_sim_in, "T", &Td);
    sim_in_set({{ ocp.model.name }}_sim_config, {{ ocp.model.name }}_sim_dims, {{ ocp.model.name }}_sim_in, "x", x_sim);
    sim_in_set({{ ocp.model.name }}_sim_config, {{ ocp.model.name }}_sim_dims, {{ ocp.model.name }}_sim_in, "u", u_sim);


    sim_status = {{ ocp.model.name }}_acados_sim_solve();
    // get and print output
    double *xn_out = calloc( {{ ocp.dims.nx }}, sizeof(double));
    sim_out_get({{ ocp.model.name }}_sim_config, {{ ocp.model.name }}_sim_dims, {{ ocp.model.name }}_sim_out, "xn", xn_out);
    printf("\nxn: \n");
    d_print_exp_mat(1, {{ ocp.dims.nx }}, xn_out, 1);

    double *S_forw_out = calloc({{ ocp.dims.nx }}*({{ ocp.dims.nx }}+{{ ocp.dims.nu }}), sizeof(double));
    if ({{ ocp.model.name }}_sim_opts->sens_forw){
        sim_out_get({{ ocp.model.name }}_sim_config, {{ ocp.model.name }}_sim_dims, {{ ocp.model.name }}_sim_out, "S_forw", S_forw_out);
        printf("\nS_forw_out: \n");
        d_print_exp_mat({{ ocp.dims.nx }}, {{ ocp.dims.nx }} + {{ ocp.dims.nu }}, S_forw_out, {{ ocp.dims.nx }});
    }

    int status = 0;
    status = acados_create();

    if (status) { 
        printf("acados_create() returned status %d. Exiting.\n", status); 
        exit(1); }

    // set initial condition
    double x0[{{ocp.dims.nx}}];
    {% for item in ocp.constraints.x0 %}
    x0[{{ loop.index0 }}] = {{ item }};
    {% endfor %}
    
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", x0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", x0);

    {% if ocp.dims.np > 0%}
    double p[{{ocp.dims.np}}];
    {% for item in ocp.constraints.p %}
    p[{{ loop.index0 }}] = {{ item }};
    {% endfor %}
    {% endif %}
    
    {% if ocp.dims.np > 0%}
    {% if ocp.solver_config.integrator_type == "IRK" %}
    for (int ii = 0; ii < {{ocp.dims.N}}; ii++) {
    impl_dae_fun[ii].set_param(impl_dae_fun+ii, p);
    impl_dae_fun_jac_x_xdot_z[ii].set_param(impl_dae_fun_jac_x_xdot_z+ii, p);
    impl_dae_jac_x_xdot_u_z[ii].set_param(impl_dae_jac_x_xdot_u_z+ii, p);
    }
    {% else %}
    for (int ii = 0; ii < {{ocp.dims.N}}; ii++) {
    forw_vde_casadi[ii].set_param(forw_vde_casadi+ii, p);
    }
    {% endif %}
    {% endif %}

    double kkt_norm_inf = 1e12, elapsed_time;

#if 1
    int NTIMINGS = 100;
    double min_time = 1e12;
    for (int ii = 0; ii < NTIMINGS; ii ++) {
        for (int i = 0; i <= nlp_dims->N; ++i) {
            blasfeo_dvecse(nlp_dims->nu[i]+nlp_dims->nx[i], 0.0, nlp_out->ux+i, 0);
        }
        status = acados_solve();
        elapsed_time = nlp_out->total_time;
        if (elapsed_time < min_time) min_time = elapsed_time;
    }
    elapsed_time = min_time;
#else
    status = acados_solve();
#endif
    kkt_norm_inf = nlp_out->inf_norm_res;
    elapsed_time = nlp_out->total_time;
    printf(" iterations %2d | time  %f |  KKT %e\n", nlp_out->sqp_iter, elapsed_time, kkt_norm_inf);

    printf("\n--- solution ---\n");
    ocp_nlp_out_print(nlp_solver->dims, nlp_out);
    if (status) { 
        printf("acados_solve() returned status %d.\n", status); 
    }

    status = acados_free();

    if (status) { 
        printf("acados_free() returned status %d. \n", status); 
    }

    return status;
}
