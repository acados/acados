/*
 * Copyright (c) The acados authors.
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

#ifndef ACADOS_SOLVER_{{ name }}_H_
#define ACADOS_SOLVER_{{ name }}_H_

#include "acados/utils/types.h"

#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"


#ifdef __cplusplus
extern "C" {
#endif

{%- if not solver_options.custom_update_filename %}
    {%- set custom_update_filename = "" %}
{% else %}
    {%- set custom_update_filename = solver_options.custom_update_filename %}
{%- endif %}

{% set cost_e = cost | last %}
{% set constraints_e = constraints | last %}
{% set dims_e = phases_dims | last %}


{% set cost_0 = cost | first %}
{% set dims_0 = phases_dims | first %}
{% set model_0 = model | first %}
{% set constraints_0 = constraints | first %}



// ** capsule for solver data **
typedef struct {{ name }}_solver_capsule
{
    // acados objects
    ocp_nlp_in *nlp_in;
    ocp_nlp_out *nlp_out;
    ocp_nlp_out *sens_out;
    ocp_nlp_solver *nlp_solver;
    void *nlp_opts;
    ocp_nlp_plan_t *nlp_solver_plan;
    ocp_nlp_config *nlp_config;
    ocp_nlp_dims *nlp_dims;

	{%- for jj in range(end=n_phases) %}{# phases loop !#}
    /* external functions phase {{ jj }} */
    // dynamics
{% if mocp_opts.integrator_type[jj] == "ERK" %}
    external_function_param_casadi *forw_vde_casadi_{{ jj }};
    external_function_param_casadi *expl_ode_fun_{{ jj }};
{% if solver_options.hessian_approx == "EXACT" %}
    external_function_param_casadi *hess_vde_casadi_{{ jj }};
{%- endif %}
{% elif mocp_opts.integrator_type[jj] == "IRK" %}
    external_function_param_{{ model[jj].dyn_ext_fun_type }} *impl_dae_fun_{{ jj }};
    external_function_param_{{ model[jj].dyn_ext_fun_type }} *impl_dae_fun_jac_x_xdot_z_{{ jj }};
    external_function_param_{{ model[jj].dyn_ext_fun_type }} *impl_dae_jac_x_xdot_u_z_{{ jj }};
{% if solver_options.hessian_approx == "EXACT" %}
    external_function_param_{{ model[jj].dyn_ext_fun_type }} *impl_dae_hess_{{ jj }};
{%- endif %}
{% elif mocp_opts.integrator_type[jj] == "LIFTED_IRK" %}
    external_function_param_{{ model[jj].dyn_ext_fun_type }} *impl_dae_fun_{{ jj }};
    external_function_param_{{ model[jj].dyn_ext_fun_type }} *impl_dae_fun_jac_x_xdot_u_{{ jj }};
{% elif mocp_opts.integrator_type[jj] == "GNSF" %}
    external_function_param_casadi *gnsf_phi_fun_{{ jj }};
    external_function_param_casadi *gnsf_phi_fun_jac_y_{{ jj }};
    external_function_param_casadi *gnsf_phi_jac_y_uhat_{{ jj }};
    external_function_param_casadi *gnsf_f_lo_jac_x1_x1dot_u_z_{{ jj }};
    external_function_param_casadi *gnsf_get_matrices_fun_{{ jj }};
{% elif mocp_opts.integrator_type[jj] == "DISCRETE" %}
    external_function_param_{{ model[jj].dyn_ext_fun_type }} *discr_dyn_phi_fun_{{ jj }};
    external_function_param_{{ model[jj].dyn_ext_fun_type }} *discr_dyn_phi_fun_jac_ut_xt_{{ jj }};
{%- if solver_options.hessian_approx == "EXACT" %}
    external_function_param_{{ model[jj].dyn_ext_fun_type }} *discr_dyn_phi_fun_jac_ut_xt_hess_{{ jj }};
{%- endif %}
{%- endif %}

    // constraints
{%- if constraints[jj].constr_type == "BGP" %}
    external_function_param_casadi *phi_constraint_{{ jj }};
{% elif constraints[jj].constr_type == "BGH" and phases_dims[jj].nh > 0 %}
    external_function_param_casadi *nl_constr_h_fun_jac_{{ jj }};
    external_function_param_casadi *nl_constr_h_fun_{{ jj }};
{%- if solver_options.hessian_approx == "EXACT" %}
    external_function_param_casadi *nl_constr_h_fun_jac_hess_{{ jj }};
{%- endif %}
{%- endif %}


    // cost
{% if cost[jj].cost_type == "NONLINEAR_LS" %}
    external_function_param_casadi *cost_y_fun_{{ jj }};
    external_function_param_casadi *cost_y_fun_jac_ut_xt_{{ jj }};
    external_function_param_casadi *cost_y_hess_{{ jj }};
{% elif cost[jj].cost_type == "CONVEX_OVER_NONLINEAR" %}
    external_function_param_casadi *conl_cost_fun_{{ jj }};
    external_function_param_casadi *conl_cost_fun_jac_hess_{{ jj }};
{%- elif cost[jj].cost_type == "EXTERNAL" %}
    external_function_param_{{ cost[jj].cost_ext_fun_type }} *ext_cost_fun_{{ jj }};
    external_function_param_{{ cost[jj].cost_ext_fun_type }} *ext_cost_fun_jac_{{ jj }};
    external_function_param_{{ cost[jj].cost_ext_fun_type }} *ext_cost_fun_jac_hess_{{ jj }};
{% endif %}
	{%- endfor %}{# for jj in range(end=n_phases) #}


{% if cost_0.cost_type_0 == "NONLINEAR_LS" %}
    external_function_param_casadi cost_y_0_fun;
    external_function_param_casadi cost_y_0_fun_jac_ut_xt;
    external_function_param_casadi cost_y_0_hess;
{% elif cost_0.cost_type_0 == "CONVEX_OVER_NONLINEAR" %}
    external_function_param_casadi conl_cost_0_fun;
    external_function_param_casadi conl_cost_0_fun_jac_hess;
{% elif cost_0.cost_type_0 == "EXTERNAL" %}
    external_function_param_{{ cost_0.cost_ext_fun_type_0 }} ext_cost_0_fun;
    external_function_param_{{ cost_0.cost_ext_fun_type_0 }} ext_cost_0_fun_jac;
    external_function_param_{{ cost_0.cost_ext_fun_type_0 }} ext_cost_0_fun_jac_hess;
{%- endif %}


{% if constraints[0].constr_type_0 == "BGP" %}
    external_function_param_casadi phi_0_constraint;
{% elif constraints[0].constr_type_0 == "BGH" and dims_0.nh_0 > 0 %}
    external_function_param_casadi nl_constr_h_0_fun_jac;
    external_function_param_casadi nl_constr_h_0_fun;
{%- if solver_options.hessian_approx == "EXACT" %}
    external_function_param_casadi nl_constr_h_0_fun_jac_hess;
{%- endif %}
{%- endif %}

{% if cost_e.cost_type_e == "NONLINEAR_LS" %}
    external_function_param_casadi cost_y_e_fun;
    external_function_param_casadi cost_y_e_fun_jac_ut_xt;
    external_function_param_casadi cost_y_e_hess;
{% elif cost_e.cost_type_e == "CONVEX_OVER_NONLINEAR" %}
    external_function_param_casadi conl_cost_e_fun;
    external_function_param_casadi conl_cost_e_fun_jac_hess;
{% elif cost_e.cost_type_e == "EXTERNAL" %}
    external_function_param_{{ cost_e.cost_ext_fun_type_e }} ext_cost_e_fun;
    external_function_param_{{ cost_e.cost_ext_fun_type_e }} ext_cost_e_fun_jac;
    external_function_param_{{ cost_e.cost_ext_fun_type_e }} ext_cost_e_fun_jac_hess;
{%- endif %}


{% if constraints_e.constr_type_e == "BGP" %}
    external_function_param_casadi phi_e_constraint;
{% elif constraints_e.constr_type_e == "BGH" and dims_e.nh_e > 0 %}
    external_function_param_casadi nl_constr_h_e_fun_jac;
    external_function_param_casadi nl_constr_h_e_fun;
{%- if solver_options.hessian_approx == "EXACT" %}
    external_function_param_casadi nl_constr_h_e_fun_jac_hess;
{%- endif %}
{%- endif %}


{%- if custom_update_filename != "" %}
void * custom_update_memory;
{%- endif %}

} {{ name }}_solver_capsule;

ACADOS_SYMBOL_EXPORT {{ name }}_solver_capsule * {{ name }}_acados_create_capsule(void);
ACADOS_SYMBOL_EXPORT int {{ name }}_acados_free_capsule({{ name }}_solver_capsule *capsule);

ACADOS_SYMBOL_EXPORT int {{ name }}_acados_create({{ name }}_solver_capsule * capsule);

ACADOS_SYMBOL_EXPORT int {{ name }}_acados_reset({{ name }}_solver_capsule* capsule, int reset_qp_solver_mem);
ACADOS_SYMBOL_EXPORT int {{ name }}_acados_create_with_discretization({{ name }}_solver_capsule* capsule, int N, double* new_time_steps);


ACADOS_SYMBOL_EXPORT int {{ name }}_acados_update_params({{ name }}_solver_capsule * capsule, int stage, double *value, int np);
ACADOS_SYMBOL_EXPORT int {{ name }}_acados_update_params_sparse({{ name }}_solver_capsule * capsule, int stage, int *idx, double *p, int n_update);

ACADOS_SYMBOL_EXPORT int {{ name }}_acados_solve({{ name }}_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT int {{ name }}_acados_free({{ name }}_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT void {{ name }}_acados_print_stats({{ name }}_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT int {{ name }}_acados_custom_update({{ name }}_solver_capsule* capsule, double* data, int data_len);


ACADOS_SYMBOL_EXPORT ocp_nlp_in *{{ name }}_acados_get_nlp_in({{ name }}_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_out *{{ name }}_acados_get_nlp_out({{ name }}_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_out *{{ name }}_acados_get_sens_out({{ name }}_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_solver *{{ name }}_acados_get_nlp_solver({{ name }}_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_config *{{ name }}_acados_get_nlp_config({{ name }}_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT void *{{ name }}_acados_get_nlp_opts({{ name }}_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_dims *{{ name }}_acados_get_nlp_dims({{ name }}_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_plan_t *{{ name }}_acados_get_nlp_plan({{ name }}_solver_capsule * capsule);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SOLVER_{{ name }}_H_


{# initial and terminal constraints and costs are not subject to multiphase #}


