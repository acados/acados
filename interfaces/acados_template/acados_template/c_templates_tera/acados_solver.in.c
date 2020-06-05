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

// example specific
#include "{{ model.name }}_model/{{ model.name }}_model.h"
{% if constraints.constr_type == "BGP" and dims.nphi %}
#include "{{ model.name }}_constraints/{{ model.name }}_phi_constraint.h"
// #include "{{ model.name }}_r_constraint/{{ model.name }}_r_constraint.h"
{% endif %}
{% if constraints.constr_type_e == "BGP" and dims.nphi_e > 0 %}
#include "{{ model.name }}_constraints/{{ model.name }}_phi_e_constraint.h"
// #include "{{ model.name }}_r_e_constraint/{{ model.name }}_r_e_constraint.h"
{% endif %}
{% if constraints.constr_type == "BGH" and dims.nh > 0 %}
#include "{{ model.name }}_constraints/{{ model.name }}_h_constraint.h"
{% endif %}
{% if constraints.constr_type_e == "BGH" and dims.nh_e > 0 %}
#include "{{ model.name }}_constraints/{{ model.name }}_h_e_constraint.h"
{% endif %}
{%- if cost.cost_type == "NONLINEAR_LS" %}
#include "{{ model.name }}_cost/{{ model.name }}_cost_y_fun.h"
{%- elif cost.cost_type == "EXTERNAL" %}
#include "{{ model.name }}_cost/{{ model.name }}_external_cost.h"
{% endif %}
{%- if cost.cost_type_e == "NONLINEAR_LS" %}
#include "{{ model.name }}_cost/{{ model.name }}_cost_y_e_fun.h"
{%- elif cost.cost_type_e == "EXTERNAL" %}
#include "{{ model.name }}_cost/{{ model.name }}_external_cost_e.h"
{% endif %}

#include "acados_solver_{{ model.name }}.h"

#define NX     {{ dims.nx }}
#define NZ     {{ dims.nz }}
#define NU     {{ dims.nu }}
#define NP     {{ dims.np }}
#define NBX    {{ dims.nbx }}
#define NBX0   {{ dims.nbx_0 }}
#define NBU    {{ dims.nbu }}
#define NSBX   {{ dims.nsbx }}
#define NSBU   {{ dims.nsbu }}
#define NSH    {{ dims.nsh }}
#define NSG    {{ dims.nsg }}
#define NSPHI  {{ dims.nsphi }}
#define NSHN   {{ dims.nsh_e }}
#define NSGN   {{ dims.nsg_e }}
#define NSPHIN {{ dims.nsphi_e }}
#define NSBXN  {{ dims.nsbx_e }}
#define NS     {{ dims.ns }}
#define NSN    {{ dims.ns_e }}
#define NG     {{ dims.ng }}
#define NBXN   {{ dims.nbx_e }}
#define NGN    {{ dims.ng_e }}
#define NY     {{ dims.ny }}
#define NYN    {{ dims.ny_e }}
#define N      {{ dims.N }}
#define NH     {{ dims.nh }}
#define NPHI   {{ dims.nphi }}
#define NHN    {{ dims.nh_e }}
#define NPHIN  {{ dims.nphi_e }}
#define NR     {{ dims.nr }}


// ** global data **
ocp_nlp_in * nlp_in;
ocp_nlp_out * nlp_out;
ocp_nlp_solver * nlp_solver;
void * nlp_opts;
ocp_nlp_plan * nlp_solver_plan;
ocp_nlp_config * nlp_config;
ocp_nlp_dims * nlp_dims;

{%- if solver_options.integrator_type == "ERK" %}
external_function_param_casadi * forw_vde_casadi;
external_function_param_casadi * expl_ode_fun;
{% if solver_options.hessian_approx == "EXACT" %}
external_function_param_casadi * hess_vde_casadi;
{%- endif %}
{%- elif solver_options.integrator_type == "IRK" %}
external_function_param_casadi * impl_dae_fun;
external_function_param_casadi * impl_dae_fun_jac_x_xdot_z;
external_function_param_casadi * impl_dae_jac_x_xdot_u_z;
{%- if solver_options.hessian_approx == "EXACT" %}
external_function_param_casadi * impl_dae_hess;
{%- endif %}
{%- elif solver_options.integrator_type == "GNSF" %}
external_function_param_casadi * gnsf_phi_fun;
external_function_param_casadi * gnsf_phi_fun_jac_y;
external_function_param_casadi * gnsf_phi_jac_y_uhat;
external_function_param_casadi * gnsf_f_lo_jac_x1_x1dot_u_z;
external_function_param_casadi * gnsf_get_matrices_fun;
{%- endif %}

{% if constraints.constr_type == "BGH" %}
external_function_param_casadi * nl_constr_h_fun;
external_function_param_casadi * nl_constr_h_fun_jac;
{%- if solver_options.hessian_approx == "EXACT" %}
external_function_param_casadi * nl_constr_h_fun_jac_hess;
{%- endif %}
{%- elif constraints.constr_type == "BGP" %}
external_function_param_casadi * phi_constraint;
// external_function_param_casadi * r_constraint;
{%- endif %}

{% if constraints.constr_type_e == "BGH" %}
external_function_param_casadi nl_constr_h_e_fun_jac;
external_function_param_casadi nl_constr_h_e_fun;
{%- if solver_options.hessian_approx == "EXACT" %}
external_function_param_casadi nl_constr_h_e_fun_jac_hess;
{%- endif %}
{%- elif constraints.constr_type_e == "BGP" %}
external_function_param_casadi phi_e_constraint;
// external_function_param_casadi r_e_constraint;
{%- endif %}

{%- if cost.cost_type == "NONLINEAR_LS" %}
external_function_param_casadi * cost_y_fun;
external_function_param_casadi * cost_y_fun_jac_ut_xt;
external_function_param_casadi * cost_y_hess;
{%- elif cost.cost_type == "EXTERNAL" %}
external_function_param_casadi * ext_cost_fun;
external_function_param_casadi * ext_cost_fun_jac;
external_function_param_casadi * ext_cost_fun_jac_hess;
{%- endif %}
{%- if cost.cost_type_e == "NONLINEAR_LS" %}
external_function_param_casadi cost_y_e_fun;
external_function_param_casadi cost_y_e_fun_jac_ut_xt;
external_function_param_casadi cost_y_e_hess;
{%- elif cost.cost_type_e == "EXTERNAL" %}
external_function_param_casadi ext_cost_e_fun;
external_function_param_casadi ext_cost_e_fun_jac;
external_function_param_casadi ext_cost_e_fun_jac_hess;
{%- endif %}


int acados_create()
{
    int status = 0;

    /************************************************
    *  plan & config
    ************************************************/
    nlp_solver_plan = ocp_nlp_plan_create(N);
    {%- if solver_options.nlp_solver_type == "SQP" %}
    nlp_solver_plan->nlp_solver = SQP;
    {% else %}
    nlp_solver_plan->nlp_solver = SQP_RTI;
    {%- endif %}

    nlp_solver_plan->ocp_qp_solver_plan.qp_solver = {{ solver_options.qp_solver }};
    for (int i = 0; i < N; i++)
        nlp_solver_plan->nlp_cost[i] = {{ cost.cost_type }};

    nlp_solver_plan->nlp_cost[N] = {{ cost.cost_type_e }};

    for (int i = 0; i < N; i++)
    {
        nlp_solver_plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
        nlp_solver_plan->sim_solver_plan[i].sim_solver = {{ solver_options.integrator_type }};
    }

    for (int i = 0; i < N; i++)
    {
        {% if constraints.constr_type == "BGP" %}
        nlp_solver_plan->nlp_constraints[i] = BGP;
        {%- else -%}
        nlp_solver_plan->nlp_constraints[i] = BGH;
        {%- endif %}
    }

    {%- if constraints.constr_type_e == "BGP" %}
    nlp_solver_plan->nlp_constraints[N] = BGP;
    {% else %}
    nlp_solver_plan->nlp_constraints[N] = BGH;
    {%- endif %}

{%- if solver_options.hessian_approx == "EXACT" %}
    {%- if solver_options.regularize_method == "NO_REGULARIZE" %}
    nlp_solver_plan->regularization = NO_REGULARIZE;
    {%- elif solver_options.regularize_method == "MIRROR" %}
    nlp_solver_plan->regularization = MIRROR;
    {%- elif solver_options.regularize_method == "PROJECT" %}
    nlp_solver_plan->regularization = PROJECT;
    {%- elif solver_options.regularize_method == "PROJECT_REDUC_HESS" %}
    nlp_solver_plan->regularization = PROJECT_REDUC_HESS;
    {%- elif solver_options.regularize_method == "CONVEXIFY" %}
    nlp_solver_plan->regularization = CONVEXIFY;
    {%- endif %}
{%- endif %}
    nlp_config = ocp_nlp_config_create(*nlp_solver_plan);


    /************************************************
    *  dimensions
    ************************************************/
    int nx[N+1];
    int nu[N+1];
    int nbx[N+1];
    int nbu[N+1];
    int nsbx[N+1];
    int nsbu[N+1];
    int nsg[N+1];
    int nsh[N+1];
    int nsphi[N+1];
    int ns[N+1];
    int ng[N+1];
    int nh[N+1];
    int nphi[N+1];
    int nz[N+1];
    int ny[N+1];
    int nr[N+1];
    int nbxe[N+1];

    for (int i = 0; i < N+1; i++)
    {
        // common
        nx[i]     = NX;
        nu[i]     = NU;
        nz[i]     = NZ;
        ns[i]     = NS;
        // cost
        ny[i]     = NY;
        // constraints
        nbx[i]    = NBX;
        nbu[i]    = NBU;
        nsbx[i]   = NSBX;
        nsbu[i]   = NSBU;
        nsg[i] = NSG;
        nsh[i]    = NSH;
        nsphi[i]  = NSPHI;
        ng[i]     = NG;
        nh[i]     = NH;
        nphi[i]   = NPHI;
        nr[i]     = NR;
        nbxe[i]   = 0;
    }

    // for initial state
    nbx[0]  = NBX0;
    nsbx[0] = 0;
    ns[0] = NS - NSBX;
    nbxe[0] = {{ dims.nbxe_0 }};

    // terminal - common
    nu[N]   = 0;
    nz[N]   = 0;
    ns[N]   = NSN;
    // cost
    ny[N]   = NYN;
    // constraint
    nbx[N]   = NBXN;
    nbu[N]   = 0;
    ng[N]    = NGN;
    nh[N]    = NHN;
    nphi[N]  = NPHIN;
    nr[N]    = {{ dims.nr_e }};

    nsbx[N]  = NSBXN;
    nsbu[N]  = 0;
    nsg[N]   = NSGN;
    nsh[N]   = NSHN;
    nsphi[N] = NSPHIN;

    /* create and set ocp_nlp_dims */
    nlp_dims = ocp_nlp_dims_create(nlp_config);

    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nx", nx);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nu", nu);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nz", nz);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "ns", ns);

    for (int i = 0; i <= N; i++)
    {
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nbx", &nbx[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nbu", &nbu[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsbx", &nsbx[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsbu", &nsbu[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "ng", &ng[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsg", &nsg[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nbxe", &nbxe[i]);
    }

    for (int i = 0; i < N; i++)
    {
        {%- if constraints.constr_type == "BGH" and dims.nh > 0 %}
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nh", &nh[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsh", &nsh[i]);
        {%- elif constraints.constr_type == "BGP" and dims.nphi > 0 %}
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nr", &nr[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nphi", &nphi[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsphi", &nsphi[i]);
        {%- endif %}
        {%- if cost.cost_type == "NONLINEAR_LS" or cost.cost_type == "LINEAR_LS" %}
        ocp_nlp_dims_set_cost(nlp_config, nlp_dims, i, "ny", &ny[i]);
        {%- endif %}
    }

    {%- if constraints.constr_type_e == "BGH" %}
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "nh", &nh[N]);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "nsh", &nsh[N]);
    {%- elif constraints.constr_type_e == "BGP" %}
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "nr", &nr[N]);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "nphi", &nphi[N]);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "nsphi", &nsphi[N]);
    {%- endif %}
    {%- if cost.cost_type_e == "NONLINEAR_LS" or cost.cost_type_e == "LINEAR_LS" %}
    ocp_nlp_dims_set_cost(nlp_config, nlp_dims, N, "ny", &ny[N]);
    {%- endif %}

{% if solver_options.integrator_type == "GNSF" -%}
    // GNSF specific dimensions
    int gnsf_nx1 = {{ dims.gnsf_nx1 }};
    int gnsf_nz1 = {{ dims.gnsf_nz1 }};
    int gnsf_nout = {{ dims.gnsf_nout }};
    int gnsf_ny = {{ dims.gnsf_ny }};
    int gnsf_nuhat = {{ dims.gnsf_nuhat }};

    for (int i = 0; i < N; i++)
    {
        if (nlp_solver_plan->sim_solver_plan[i].sim_solver == GNSF)
        {
            ocp_nlp_dims_set_dynamics(nlp_config, nlp_dims, i, "gnsf_nx1", &gnsf_nx1);
            ocp_nlp_dims_set_dynamics(nlp_config, nlp_dims, i, "gnsf_nz1", &gnsf_nz1);
            ocp_nlp_dims_set_dynamics(nlp_config, nlp_dims, i, "gnsf_nout", &gnsf_nout);
            ocp_nlp_dims_set_dynamics(nlp_config, nlp_dims, i, "gnsf_ny", &gnsf_ny);
            ocp_nlp_dims_set_dynamics(nlp_config, nlp_dims, i, "gnsf_nuhat", &gnsf_nuhat);
        }
    }
{%- endif %}

    /************************************************
    *  external functions
    ************************************************/
    {%- if constraints.constr_type == "BGP" %}
    phi_constraint = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++)
    {
        // nonlinear part of convex-composite constraint
        phi_constraint[i].casadi_fun = &{{ model.name }}_phi_constraint;
        phi_constraint[i].casadi_n_in = &{{ model.name }}_phi_constraint_n_in;
        phi_constraint[i].casadi_n_out = &{{ model.name }}_phi_constraint_n_out;
        phi_constraint[i].casadi_sparsity_in = &{{ model.name }}_phi_constraint_sparsity_in;
        phi_constraint[i].casadi_sparsity_out = &{{ model.name }}_phi_constraint_sparsity_out;
        phi_constraint[i].casadi_work = &{{ model.name }}_phi_constraint_work;

        external_function_param_casadi_create(&phi_constraint[i], {{ dims.np }});
    }
    // r_constraint = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    // for (int i = 0; i < N; i++) {
    //     // nonlinear part of convex-composite constraint
    //     r_constraint[i].casadi_fun = &{{ model.name }}_r_constraint;
    //     r_constraint[i].casadi_n_in = &{{ model.name }}_r_constraint_n_in;
    //     r_constraint[i].casadi_n_out = &{{ model.name }}_r_constraint_n_out;
    //     r_constraint[i].casadi_sparsity_in = &{{ model.name }}_r_constraint_sparsity_in;
    //     r_constraint[i].casadi_sparsity_out = &{{ model.name }}_r_constraint_sparsity_out;
    //     r_constraint[i].casadi_work = &{{ model.name }}_r_constraint_work;

    //     external_function_param_casadi_create(&r_constraint[i], {{ dims.np }});
    // }
    {%- endif %}

    {%- if constraints.constr_type_e == "BGP" %}
    // nonlinear part of convex-composite constraint
    phi_e_constraint.casadi_fun = &{{ model.name }}_phi_e_constraint;
    phi_e_constraint.casadi_n_in = &{{ model.name }}_phi_e_constraint_n_in;
    phi_e_constraint.casadi_n_out = &{{ model.name }}_phi_e_constraint_n_out;
    phi_e_constraint.casadi_sparsity_in = &{{ model.name }}_phi_e_constraint_sparsity_in;
    phi_e_constraint.casadi_sparsity_out = &{{ model.name }}_phi_e_constraint_sparsity_out;
    phi_e_constraint.casadi_work = &{{ model.name }}_phi_e_constraint_work;

    external_function_param_casadi_create(&phi_e_constraint, {{ dims.np }});
    
    // nonlinear part of convex-composite constraint
    // r_e_constraint.casadi_fun = &{{ model.name }}_r_e_constraint;
    // r_e_constraint.casadi_n_in = &{{ model.name }}_r_e_constraint_n_in;
    // r_e_constraint.casadi_n_out = &{{ model.name }}_r_e_constraint_n_out;
    // r_e_constraint.casadi_sparsity_in = &{{ model.name }}_r_e_constraint_sparsity_in;
    // r_e_constraint.casadi_sparsity_out = &{{ model.name }}_r_e_constraint_sparsity_out;
    // r_e_constraint.casadi_work = &{{ model.name }}_r_e_constraint_work;

    // external_function_param_casadi_create(&r_e_constraint, {{ dims.np }});
    {% endif %}

    {%- if constraints.constr_type == "BGH" and dims.nh > 0  %}
    nl_constr_h_fun_jac = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        nl_constr_h_fun_jac[i].casadi_fun = &{{ model.name }}_constr_h_fun_jac_uxt_zt;
        nl_constr_h_fun_jac[i].casadi_n_in = &{{ model.name }}_constr_h_fun_jac_uxt_zt_n_in;
        nl_constr_h_fun_jac[i].casadi_n_out = &{{ model.name }}_constr_h_fun_jac_uxt_zt_n_out;
        nl_constr_h_fun_jac[i].casadi_sparsity_in = &{{ model.name }}_constr_h_fun_jac_uxt_zt_sparsity_in;
        nl_constr_h_fun_jac[i].casadi_sparsity_out = &{{ model.name }}_constr_h_fun_jac_uxt_zt_sparsity_out;
        nl_constr_h_fun_jac[i].casadi_work = &{{ model.name }}_constr_h_fun_jac_uxt_zt_work;
        external_function_param_casadi_create(&nl_constr_h_fun_jac[i], {{ dims.np }});
    }
    nl_constr_h_fun = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        nl_constr_h_fun[i].casadi_fun = &{{ model.name }}_constr_h_fun;
        nl_constr_h_fun[i].casadi_n_in = &{{ model.name }}_constr_h_fun_n_in;
        nl_constr_h_fun[i].casadi_n_out = &{{ model.name }}_constr_h_fun_n_out;
        nl_constr_h_fun[i].casadi_sparsity_in = &{{ model.name }}_constr_h_fun_sparsity_in;
        nl_constr_h_fun[i].casadi_sparsity_out = &{{ model.name }}_constr_h_fun_sparsity_out;
        nl_constr_h_fun[i].casadi_work = &{{ model.name }}_constr_h_fun_work;
        external_function_param_casadi_create(&nl_constr_h_fun[i], {{ dims.np }});
    }
    {% if solver_options.hessian_approx == "EXACT" %}
    nl_constr_h_fun_jac_hess = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        // nonlinear constraint
        nl_constr_h_fun_jac_hess[i].casadi_fun = &{{ model.name }}_constr_h_fun_jac_uxt_hess;
        nl_constr_h_fun_jac_hess[i].casadi_n_in = &{{ model.name }}_constr_h_fun_jac_uxt_hess_n_in;
        nl_constr_h_fun_jac_hess[i].casadi_n_out = &{{ model.name }}_constr_h_fun_jac_uxt_hess_n_out;
        nl_constr_h_fun_jac_hess[i].casadi_sparsity_in = &{{ model.name }}_constr_h_fun_jac_uxt_hess_sparsity_in;
        nl_constr_h_fun_jac_hess[i].casadi_sparsity_out = &{{ model.name }}_constr_h_fun_jac_uxt_hess_sparsity_out;
        nl_constr_h_fun_jac_hess[i].casadi_work = &{{ model.name }}_constr_h_fun_jac_uxt_hess_work;

        external_function_param_casadi_create(&nl_constr_h_fun_jac_hess[i], {{ dims.np }});
    }
    {% endif %}
    {% endif %}

    {%- if constraints.constr_type_e == "BGH" and dims.nh_e > 0 %}
    nl_constr_h_e_fun_jac.casadi_fun = &{{ model.name }}_constr_h_e_fun_jac_uxt_zt;
    nl_constr_h_e_fun_jac.casadi_n_in = &{{ model.name }}_constr_h_e_fun_jac_uxt_zt_n_in;
    nl_constr_h_e_fun_jac.casadi_n_out = &{{ model.name }}_constr_h_e_fun_jac_uxt_zt_n_out;
    nl_constr_h_e_fun_jac.casadi_sparsity_in = &{{ model.name }}_constr_h_e_fun_jac_uxt_zt_sparsity_in;
    nl_constr_h_e_fun_jac.casadi_sparsity_out = &{{ model.name }}_constr_h_e_fun_jac_uxt_zt_sparsity_out;
    nl_constr_h_e_fun_jac.casadi_work = &{{ model.name }}_constr_h_e_fun_jac_uxt_zt_work;
    external_function_param_casadi_create(&nl_constr_h_e_fun_jac, {{ dims.np }});

    nl_constr_h_e_fun.casadi_fun = &{{ model.name }}_constr_h_e_fun;
    nl_constr_h_e_fun.casadi_n_in = &{{ model.name }}_constr_h_e_fun_n_in;
    nl_constr_h_e_fun.casadi_n_out = &{{ model.name }}_constr_h_e_fun_n_out;
    nl_constr_h_e_fun.casadi_sparsity_in = &{{ model.name }}_constr_h_e_fun_sparsity_in;
    nl_constr_h_e_fun.casadi_sparsity_out = &{{ model.name }}_constr_h_e_fun_sparsity_out;
    nl_constr_h_e_fun.casadi_work = &{{ model.name }}_constr_h_e_fun_work;
    external_function_param_casadi_create(&nl_constr_h_e_fun, {{ dims.np }});

    {% if solver_options.hessian_approx == "EXACT" %}
    // nonlinear constraint
    nl_constr_h_e_fun_jac_hess.casadi_fun = &{{ model.name }}_constr_h_e_fun_jac_uxt_hess;
    nl_constr_h_e_fun_jac_hess.casadi_n_in = &{{ model.name }}_constr_h_e_fun_jac_uxt_hess_n_in;
    nl_constr_h_e_fun_jac_hess.casadi_n_out = &{{ model.name }}_constr_h_e_fun_jac_uxt_hess_n_out;
    nl_constr_h_e_fun_jac_hess.casadi_sparsity_in = &{{ model.name }}_constr_h_e_fun_jac_uxt_hess_sparsity_in;
    nl_constr_h_e_fun_jac_hess.casadi_sparsity_out = &{{ model.name }}_constr_h_e_fun_jac_uxt_hess_sparsity_out;
    nl_constr_h_e_fun_jac_hess.casadi_work = &{{ model.name }}_constr_h_e_fun_jac_uxt_hess_work;
    external_function_param_casadi_create(&nl_constr_h_e_fun_jac_hess, {{ dims.np }});
    {% endif %}
    {%- endif %}

{% if solver_options.integrator_type == "ERK" %}
    // explicit ode
    forw_vde_casadi = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        forw_vde_casadi[i].casadi_fun = &{{ model.name }}_expl_vde_forw;
        forw_vde_casadi[i].casadi_n_in = &{{ model.name }}_expl_vde_forw_n_in;
        forw_vde_casadi[i].casadi_n_out = &{{ model.name }}_expl_vde_forw_n_out;
        forw_vde_casadi[i].casadi_sparsity_in = &{{ model.name }}_expl_vde_forw_sparsity_in;
        forw_vde_casadi[i].casadi_sparsity_out = &{{ model.name }}_expl_vde_forw_sparsity_out;
        forw_vde_casadi[i].casadi_work = &{{ model.name }}_expl_vde_forw_work;
        external_function_param_casadi_create(&forw_vde_casadi[i], {{ dims.np }});
    }

    expl_ode_fun = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        expl_ode_fun[i].casadi_fun = &{{ model.name }}_expl_ode_fun;
        expl_ode_fun[i].casadi_n_in = &{{ model.name }}_expl_ode_fun_n_in;
        expl_ode_fun[i].casadi_n_out = &{{ model.name }}_expl_ode_fun_n_out;
        expl_ode_fun[i].casadi_sparsity_in = &{{ model.name }}_expl_ode_fun_sparsity_in;
        expl_ode_fun[i].casadi_sparsity_out = &{{ model.name }}_expl_ode_fun_sparsity_out;
        expl_ode_fun[i].casadi_work = &{{ model.name }}_expl_ode_fun_work;
        external_function_param_casadi_create(&expl_ode_fun[i], {{ dims.np }});
    }

    {%- if solver_options.hessian_approx == "EXACT" %}
    hess_vde_casadi = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        hess_vde_casadi[i].casadi_fun = &{{ model.name }}_expl_ode_hess;
        hess_vde_casadi[i].casadi_n_in = &{{ model.name }}_expl_ode_hess_n_in;
        hess_vde_casadi[i].casadi_n_out = &{{ model.name }}_expl_ode_hess_n_out;
        hess_vde_casadi[i].casadi_sparsity_in = &{{ model.name }}_expl_ode_hess_sparsity_in;
        hess_vde_casadi[i].casadi_sparsity_out = &{{ model.name }}_expl_ode_hess_sparsity_out;
        hess_vde_casadi[i].casadi_work = &{{ model.name }}_expl_ode_hess_work;
        external_function_param_casadi_create(&hess_vde_casadi[i], {{ dims.np }});
    }
    {%- endif %}

{% elif solver_options.integrator_type == "IRK" %}
    // implicit dae
    impl_dae_fun = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        impl_dae_fun[i].casadi_fun = &{{ model.name }}_impl_dae_fun;
        impl_dae_fun[i].casadi_work = &{{ model.name }}_impl_dae_fun_work;
        impl_dae_fun[i].casadi_sparsity_in = &{{ model.name }}_impl_dae_fun_sparsity_in;
        impl_dae_fun[i].casadi_sparsity_out = &{{ model.name }}_impl_dae_fun_sparsity_out;
        impl_dae_fun[i].casadi_n_in = &{{ model.name }}_impl_dae_fun_n_in;
        impl_dae_fun[i].casadi_n_out = &{{ model.name }}_impl_dae_fun_n_out;
        external_function_param_casadi_create(&impl_dae_fun[i], {{ dims.np }});
    }

    impl_dae_fun_jac_x_xdot_z = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        impl_dae_fun_jac_x_xdot_z[i].casadi_fun = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z;
        impl_dae_fun_jac_x_xdot_z[i].casadi_work = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z_work;
        impl_dae_fun_jac_x_xdot_z[i].casadi_sparsity_in = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z_sparsity_in;
        impl_dae_fun_jac_x_xdot_z[i].casadi_sparsity_out = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z_sparsity_out;
        impl_dae_fun_jac_x_xdot_z[i].casadi_n_in = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z_n_in;
        impl_dae_fun_jac_x_xdot_z[i].casadi_n_out = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z_n_out;
        external_function_param_casadi_create(&impl_dae_fun_jac_x_xdot_z[i], {{ dims.np }});
    }

    impl_dae_jac_x_xdot_u_z = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        impl_dae_jac_x_xdot_u_z[i].casadi_fun = &{{ model.name }}_impl_dae_jac_x_xdot_u_z;
        impl_dae_jac_x_xdot_u_z[i].casadi_work = &{{ model.name }}_impl_dae_jac_x_xdot_u_z_work;
        impl_dae_jac_x_xdot_u_z[i].casadi_sparsity_in = &{{ model.name }}_impl_dae_jac_x_xdot_u_z_sparsity_in;
        impl_dae_jac_x_xdot_u_z[i].casadi_sparsity_out = &{{ model.name }}_impl_dae_jac_x_xdot_u_z_sparsity_out;
        impl_dae_jac_x_xdot_u_z[i].casadi_n_in = &{{ model.name }}_impl_dae_jac_x_xdot_u_z_n_in;
        impl_dae_jac_x_xdot_u_z[i].casadi_n_out = &{{ model.name }}_impl_dae_jac_x_xdot_u_z_n_out;
        external_function_param_casadi_create(&impl_dae_jac_x_xdot_u_z[i], {{ dims.np }});
    }

    {%- if solver_options.hessian_approx == "EXACT" %}
    impl_dae_hess = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        impl_dae_hess[i].casadi_fun = &{{ model.name }}_impl_dae_hess;
        impl_dae_hess[i].casadi_work = &{{ model.name }}_impl_dae_hess_work;
        impl_dae_hess[i].casadi_sparsity_in = &{{ model.name }}_impl_dae_hess_sparsity_in;
        impl_dae_hess[i].casadi_sparsity_out = &{{ model.name }}_impl_dae_hess_sparsity_out;
        impl_dae_hess[i].casadi_n_in = &{{ model.name }}_impl_dae_hess_n_in;
        impl_dae_hess[i].casadi_n_out = &{{ model.name }}_impl_dae_hess_n_out;
        external_function_param_casadi_create(&impl_dae_hess[i], {{ dims.np }});
    }
    {%- endif %}

{% elif solver_options.integrator_type == "GNSF" %}
    gnsf_phi_fun = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        gnsf_phi_fun[i].casadi_fun = &{{ model.name }}_gnsf_phi_fun;
        gnsf_phi_fun[i].casadi_work = &{{ model.name }}_gnsf_phi_fun_work;
        gnsf_phi_fun[i].casadi_sparsity_in = &{{ model.name }}_gnsf_phi_fun_sparsity_in;
        gnsf_phi_fun[i].casadi_sparsity_out = &{{ model.name }}_gnsf_phi_fun_sparsity_out;
        gnsf_phi_fun[i].casadi_n_in = &{{ model.name }}_gnsf_phi_fun_n_in;
        gnsf_phi_fun[i].casadi_n_out = &{{ model.name }}_gnsf_phi_fun_n_out;
        external_function_param_casadi_create(&gnsf_phi_fun[i], {{ dims.np }});
    }

    gnsf_phi_fun_jac_y = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        gnsf_phi_fun_jac_y[i].casadi_fun = &{{ model.name }}_gnsf_phi_fun_jac_y;
        gnsf_phi_fun_jac_y[i].casadi_work = &{{ model.name }}_gnsf_phi_fun_jac_y_work;
        gnsf_phi_fun_jac_y[i].casadi_sparsity_in = &{{ model.name }}_gnsf_phi_fun_jac_y_sparsity_in;
        gnsf_phi_fun_jac_y[i].casadi_sparsity_out = &{{ model.name }}_gnsf_phi_fun_jac_y_sparsity_out;
        gnsf_phi_fun_jac_y[i].casadi_n_in = &{{ model.name }}_gnsf_phi_fun_jac_y_n_in;
        gnsf_phi_fun_jac_y[i].casadi_n_out = &{{ model.name }}_gnsf_phi_fun_jac_y_n_out;
        external_function_param_casadi_create(&gnsf_phi_fun_jac_y[i], {{ dims.np }});
    }

    gnsf_phi_jac_y_uhat = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        gnsf_phi_jac_y_uhat[i].casadi_fun = &{{ model.name }}_gnsf_phi_jac_y_uhat;
        gnsf_phi_jac_y_uhat[i].casadi_work = &{{ model.name }}_gnsf_phi_jac_y_uhat_work;
        gnsf_phi_jac_y_uhat[i].casadi_sparsity_in = &{{ model.name }}_gnsf_phi_jac_y_uhat_sparsity_in;
        gnsf_phi_jac_y_uhat[i].casadi_sparsity_out = &{{ model.name }}_gnsf_phi_jac_y_uhat_sparsity_out;
        gnsf_phi_jac_y_uhat[i].casadi_n_in = &{{ model.name }}_gnsf_phi_jac_y_uhat_n_in;
        gnsf_phi_jac_y_uhat[i].casadi_n_out = &{{ model.name }}_gnsf_phi_jac_y_uhat_n_out;
        external_function_param_casadi_create(&gnsf_phi_jac_y_uhat[i], {{ dims.np }});
    }

    gnsf_f_lo_jac_x1_x1dot_u_z = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        gnsf_f_lo_jac_x1_x1dot_u_z[i].casadi_fun = &{{ model.name }}_gnsf_f_lo_fun_jac_x1k1uz;
        gnsf_f_lo_jac_x1_x1dot_u_z[i].casadi_work = &{{ model.name }}_gnsf_f_lo_fun_jac_x1k1uz_work;
        gnsf_f_lo_jac_x1_x1dot_u_z[i].casadi_sparsity_in = &{{ model.name }}_gnsf_f_lo_fun_jac_x1k1uz_sparsity_in;
        gnsf_f_lo_jac_x1_x1dot_u_z[i].casadi_sparsity_out = &{{ model.name }}_gnsf_f_lo_fun_jac_x1k1uz_sparsity_out;
        gnsf_f_lo_jac_x1_x1dot_u_z[i].casadi_n_in = &{{ model.name }}_gnsf_f_lo_fun_jac_x1k1uz_n_in;
        gnsf_f_lo_jac_x1_x1dot_u_z[i].casadi_n_out = &{{ model.name }}_gnsf_f_lo_fun_jac_x1k1uz_n_out;
        external_function_param_casadi_create(&gnsf_f_lo_jac_x1_x1dot_u_z[i], {{ dims.np }});
    }

    gnsf_get_matrices_fun = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++) {
        gnsf_get_matrices_fun[i].casadi_fun = &{{ model.name }}_gnsf_get_matrices_fun;
        gnsf_get_matrices_fun[i].casadi_work = &{{ model.name }}_gnsf_get_matrices_fun_work;
        gnsf_get_matrices_fun[i].casadi_sparsity_in = &{{ model.name }}_gnsf_get_matrices_fun_sparsity_in;
        gnsf_get_matrices_fun[i].casadi_sparsity_out = &{{ model.name }}_gnsf_get_matrices_fun_sparsity_out;
        gnsf_get_matrices_fun[i].casadi_n_in = &{{ model.name }}_gnsf_get_matrices_fun_n_in;
        gnsf_get_matrices_fun[i].casadi_n_out = &{{ model.name }}_gnsf_get_matrices_fun_n_out;
        external_function_param_casadi_create(&gnsf_get_matrices_fun[i], {{ dims.np }});
    }
{%- endif %}

{%- if cost.cost_type == "NONLINEAR_LS" %}
    // nonlinear least squares cost
    cost_y_fun = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++)
    {
        cost_y_fun[i].casadi_fun = &{{ model.name }}_cost_y_fun;
        cost_y_fun[i].casadi_n_in = &{{ model.name }}_cost_y_fun_n_in;
        cost_y_fun[i].casadi_n_out = &{{ model.name }}_cost_y_fun_n_out;
        cost_y_fun[i].casadi_sparsity_in = &{{ model.name }}_cost_y_fun_sparsity_in;
        cost_y_fun[i].casadi_sparsity_out = &{{ model.name }}_cost_y_fun_sparsity_out;
        cost_y_fun[i].casadi_work = &{{ model.name }}_cost_y_fun_work;

        external_function_param_casadi_create(&cost_y_fun[i], {{ dims.np }});
    }

    cost_y_fun_jac_ut_xt = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++)
    {
        cost_y_fun_jac_ut_xt[i].casadi_fun = &{{ model.name }}_cost_y_fun_jac_ut_xt;
        cost_y_fun_jac_ut_xt[i].casadi_n_in = &{{ model.name }}_cost_y_fun_jac_ut_xt_n_in;
        cost_y_fun_jac_ut_xt[i].casadi_n_out = &{{ model.name }}_cost_y_fun_jac_ut_xt_n_out;
        cost_y_fun_jac_ut_xt[i].casadi_sparsity_in = &{{ model.name }}_cost_y_fun_jac_ut_xt_sparsity_in;
        cost_y_fun_jac_ut_xt[i].casadi_sparsity_out = &{{ model.name }}_cost_y_fun_jac_ut_xt_sparsity_out;
        cost_y_fun_jac_ut_xt[i].casadi_work = &{{ model.name }}_cost_y_fun_jac_ut_xt_work;

        external_function_param_casadi_create(&cost_y_fun_jac_ut_xt[i], {{ dims.np }});
    }

    cost_y_hess = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++)
    {
        cost_y_hess[i].casadi_fun = &{{ model.name }}_cost_y_hess;
        cost_y_hess[i].casadi_n_in = &{{ model.name }}_cost_y_hess_n_in;
        cost_y_hess[i].casadi_n_out = &{{ model.name }}_cost_y_hess_n_out;
        cost_y_hess[i].casadi_sparsity_in = &{{ model.name }}_cost_y_hess_sparsity_in;
        cost_y_hess[i].casadi_sparsity_out = &{{ model.name }}_cost_y_hess_sparsity_out;
        cost_y_hess[i].casadi_work = &{{ model.name }}_cost_y_hess_work;

        external_function_param_casadi_create(&cost_y_hess[i], {{ dims.np }});
    }
{%- elif cost.cost_type == "EXTERNAL" %}
    // external cost
    ext_cost_fun = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++)
    {
        ext_cost_fun[i].casadi_fun = &{{ model.name }}_ext_cost_fun;
        ext_cost_fun[i].casadi_n_in = &{{ model.name }}_ext_cost_fun_n_in;
        ext_cost_fun[i].casadi_n_out = &{{ model.name }}_ext_cost_fun_n_out;
        ext_cost_fun[i].casadi_sparsity_in = &{{ model.name }}_ext_cost_fun_sparsity_in;
        ext_cost_fun[i].casadi_sparsity_out = &{{ model.name }}_ext_cost_fun_sparsity_out;
        ext_cost_fun[i].casadi_work = &{{ model.name }}_ext_cost_fun_work;

        external_function_param_casadi_create(&ext_cost_fun[i], {{ dims.np }});
    }

    ext_cost_fun_jac = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++)
    {
        // residual function
        ext_cost_fun_jac[i].casadi_fun = &{{ model.name }}_ext_cost_fun_jac;
        ext_cost_fun_jac[i].casadi_n_in = &{{ model.name }}_ext_cost_fun_jac_n_in;
        ext_cost_fun_jac[i].casadi_n_out = &{{ model.name }}_ext_cost_fun_jac_n_out;
        ext_cost_fun_jac[i].casadi_sparsity_in = &{{ model.name }}_ext_cost_fun_jac_sparsity_in;
        ext_cost_fun_jac[i].casadi_sparsity_out = &{{ model.name }}_ext_cost_fun_jac_sparsity_out;
        ext_cost_fun_jac[i].casadi_work = &{{ model.name }}_ext_cost_fun_jac_work;

        external_function_param_casadi_create(&ext_cost_fun_jac[i], {{ dims.np }});
    }

    ext_cost_fun_jac_hess = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; i++)
    {
        // residual function
        ext_cost_fun_jac_hess[i].casadi_fun = &{{ model.name }}_ext_cost_fun_jac_hess;
        ext_cost_fun_jac_hess[i].casadi_n_in = &{{ model.name }}_ext_cost_fun_jac_hess_n_in;
        ext_cost_fun_jac_hess[i].casadi_n_out = &{{ model.name }}_ext_cost_fun_jac_hess_n_out;
        ext_cost_fun_jac_hess[i].casadi_sparsity_in = &{{ model.name }}_ext_cost_fun_jac_hess_sparsity_in;
        ext_cost_fun_jac_hess[i].casadi_sparsity_out = &{{ model.name }}_ext_cost_fun_jac_hess_sparsity_out;
        ext_cost_fun_jac_hess[i].casadi_work = &{{ model.name }}_ext_cost_fun_jac_hess_work;

        external_function_param_casadi_create(&ext_cost_fun_jac_hess[i], {{ dims.np }});
    }
{%- endif %}

{%- if cost.cost_type_e == "NONLINEAR_LS" %}
    // nonlinear least square function
    cost_y_e_fun.casadi_fun = &{{ model.name }}_cost_y_e_fun;
    cost_y_e_fun.casadi_n_in = &{{ model.name }}_cost_y_e_fun_n_in;
    cost_y_e_fun.casadi_n_out = &{{ model.name }}_cost_y_e_fun_n_out;
    cost_y_e_fun.casadi_sparsity_in = &{{ model.name }}_cost_y_e_fun_sparsity_in;
    cost_y_e_fun.casadi_sparsity_out = &{{ model.name }}_cost_y_e_fun_sparsity_out;
    cost_y_e_fun.casadi_work = &{{ model.name }}_cost_y_e_fun_work;
    external_function_param_casadi_create(&cost_y_e_fun, {{ dims.np }});

    cost_y_e_fun_jac_ut_xt.casadi_fun = &{{ model.name }}_cost_y_e_fun_jac_ut_xt;
    cost_y_e_fun_jac_ut_xt.casadi_n_in = &{{ model.name }}_cost_y_e_fun_jac_ut_xt_n_in;
    cost_y_e_fun_jac_ut_xt.casadi_n_out = &{{ model.name }}_cost_y_e_fun_jac_ut_xt_n_out;
    cost_y_e_fun_jac_ut_xt.casadi_sparsity_in = &{{ model.name }}_cost_y_e_fun_jac_ut_xt_sparsity_in;
    cost_y_e_fun_jac_ut_xt.casadi_sparsity_out = &{{ model.name }}_cost_y_e_fun_jac_ut_xt_sparsity_out;
    cost_y_e_fun_jac_ut_xt.casadi_work = &{{ model.name }}_cost_y_e_fun_jac_ut_xt_work;
    external_function_param_casadi_create(&cost_y_e_fun_jac_ut_xt, {{ dims.np }});

    cost_y_e_hess.casadi_fun = &{{ model.name }}_cost_y_e_hess;
    cost_y_e_hess.casadi_n_in = &{{ model.name }}_cost_y_e_hess_n_in;
    cost_y_e_hess.casadi_n_out = &{{ model.name }}_cost_y_e_hess_n_out;
    cost_y_e_hess.casadi_sparsity_in = &{{ model.name }}_cost_y_e_hess_sparsity_in;
    cost_y_e_hess.casadi_sparsity_out = &{{ model.name }}_cost_y_e_hess_sparsity_out;
    cost_y_e_hess.casadi_work = &{{ model.name }}_cost_y_e_hess_work;
    external_function_param_casadi_create(&cost_y_e_hess, {{ dims.np }});

{%- elif cost.cost_type_e == "EXTERNAL" %}
    // external cost
    ext_cost_e_fun.casadi_fun = &{{ model.name }}_ext_cost_e_fun;
    ext_cost_e_fun.casadi_n_in = &{{ model.name }}_ext_cost_e_fun_n_in;
    ext_cost_e_fun.casadi_n_out = &{{ model.name }}_ext_cost_e_fun_n_out;
    ext_cost_e_fun.casadi_sparsity_in = &{{ model.name }}_ext_cost_e_fun_sparsity_in;
    ext_cost_e_fun.casadi_sparsity_out = &{{ model.name }}_ext_cost_e_fun_sparsity_out;
    ext_cost_e_fun.casadi_work = &{{ model.name }}_ext_cost_e_fun_work;
    external_function_param_casadi_create(&ext_cost_e_fun, {{ dims.np }});

    // external cost
    ext_cost_e_fun_jac.casadi_fun = &{{ model.name }}_ext_cost_e_fun_jac;
    ext_cost_e_fun_jac.casadi_n_in = &{{ model.name }}_ext_cost_e_fun_jac_n_in;
    ext_cost_e_fun_jac.casadi_n_out = &{{ model.name }}_ext_cost_e_fun_jac_n_out;
    ext_cost_e_fun_jac.casadi_sparsity_in = &{{ model.name }}_ext_cost_e_fun_jac_sparsity_in;
    ext_cost_e_fun_jac.casadi_sparsity_out = &{{ model.name }}_ext_cost_e_fun_jac_sparsity_out;
    ext_cost_e_fun_jac.casadi_work = &{{ model.name }}_ext_cost_e_fun_jac_work;
    external_function_param_casadi_create(&ext_cost_e_fun_jac, {{ dims.np }});

    // external cost
    ext_cost_e_fun_jac_hess.casadi_fun = &{{ model.name }}_ext_cost_e_fun_jac_hess;
    ext_cost_e_fun_jac_hess.casadi_n_in = &{{ model.name }}_ext_cost_e_fun_jac_hess_n_in;
    ext_cost_e_fun_jac_hess.casadi_n_out = &{{ model.name }}_ext_cost_e_fun_jac_hess_n_out;
    ext_cost_e_fun_jac_hess.casadi_sparsity_in = &{{ model.name }}_ext_cost_e_fun_jac_hess_sparsity_in;
    ext_cost_e_fun_jac_hess.casadi_sparsity_out = &{{ model.name }}_ext_cost_e_fun_jac_hess_sparsity_out;
    ext_cost_e_fun_jac_hess.casadi_work = &{{ model.name }}_ext_cost_e_fun_jac_hess_work;
    external_function_param_casadi_create(&ext_cost_e_fun_jac_hess, {{ dims.np }});
{%- endif %}

    /************************************************
    *  nlp_in
    ************************************************/
    nlp_in = ocp_nlp_in_create(nlp_config, nlp_dims);

    double time_steps[N];
    {%- for j in range(end=dims.N) %}
    time_steps[{{ j }}] = {{ solver_options.time_steps[j] }};
    {%- endfor %}

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_in_set(nlp_config, nlp_dims, nlp_in, i, "Ts", &time_steps[i]);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "scaling", &time_steps[i]);
    }

    /**** Dynamics ****/
    for (int i = 0; i < N; i++)
    {
    {%- if solver_options.integrator_type == "ERK" %}
        ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "expl_vde_forw", &forw_vde_casadi[i]);
        ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "expl_ode_fun", &expl_ode_fun[i]);
        {%- if solver_options.hessian_approx == "EXACT" %}
        ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "expl_ode_hess", &hess_vde_casadi[i]);
        {%- endif %}
    {% elif solver_options.integrator_type == "IRK" %}
        ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "impl_dae_fun", &impl_dae_fun[i]);
        ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i,
                                   "impl_dae_fun_jac_x_xdot_z", &impl_dae_fun_jac_x_xdot_z[i]);
        ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i,
                                   "impl_dae_jac_x_xdot_u", &impl_dae_jac_x_xdot_u_z[i]);
        {%- if solver_options.hessian_approx == "EXACT" %}
        ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "impl_dae_hess", &impl_dae_hess[i]);
        {%- endif %}
    {% elif solver_options.integrator_type == "GNSF" %}
        ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "phi_fun", &gnsf_phi_fun[i]);
        ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "phi_fun_jac_y", &gnsf_phi_fun_jac_y[i]);
        ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "phi_jac_y_uhat", &gnsf_phi_jac_y_uhat[i]);
        ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "f_lo_jac_x1_x1dot_u_z",
                                   &gnsf_f_lo_jac_x1_x1dot_u_z[i]);
        ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "gnsf_get_matrices_fun",
                                   &gnsf_get_matrices_fun[i]);
    {%- endif %}
    }


    /**** Cost ****/
{%- if cost.cost_type == "NONLINEAR_LS" or cost.cost_type == "LINEAR_LS" %}
{% if dims.ny > 0 %}
    double W[NY*NY];
    {% for j in range(end=dims.ny) %}
        {%- for k in range(end=dims.ny) %}
    W[{{ j }}+(NY) * {{ k }}] = {{ cost.W[j][k] }};
        {%- endfor %}
    {%- endfor %}

    double yref[NY];
    {% for j in range(end=dims.ny) %}
    yref[{{ j }}] = {{ cost.yref[j] }};
    {%- endfor %}

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "W", W);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "yref", yref);
    }
{% endif %}
{% endif %}

{%- if cost.cost_type == "LINEAR_LS" %}
    double Vx[NY*NX];
    {% for j in range(end=dims.ny) %}
        {%- for k in range(end=dims.nx) %}
    Vx[{{ j }}+(NY) * {{ k }}] = {{ cost.Vx[j][k] }};
        {%- endfor %}
    {%- endfor %}
    for (int i = 0; i < N; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Vx", Vx);
    }

{% if dims.ny > 0 and dims.nu > 0 %}
    double Vu[NY*NU];
    {% for j in range(end=dims.ny) %}
        {%- for k in range(end=dims.nu) %}
    Vu[{{ j }}+(NY) * {{ k }}] = {{ cost.Vu[j][k] }};
        {%- endfor %}
    {%- endfor %}

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Vu", Vu);
    }
{% endif %}

{% if dims.ny > 0 and dims.nz > 0 %}
    double Vz[NY*NZ];
    {% for j in range(end=dims.ny) %}
        {%- for k in range(end=dims.nz) %}
    Vz[{{ j }}+(NY) * {{ k }}] = {{ cost.Vz[j][k] }};
        {%- endfor %}
    {%- endfor %}

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Vz", Vz);
    }
{%- endif %}
{%- endif %}{# LINEAR LS #}

{%- if cost.cost_type == "NONLINEAR_LS" %}
    for (int i = 0; i < N; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun", &cost_y_fun[i]);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun_jac", &cost_y_fun_jac_ut_xt[i]);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "nls_y_hess", &cost_y_hess[i]);
    }
{%- elif cost.cost_type == "EXTERNAL" %}
    for (int i = 0; i < N; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "ext_cost_fun", &ext_cost_fun[i]);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "ext_cost_fun_jac", &ext_cost_fun_jac[i]);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "ext_cost_fun_jac_hess", &ext_cost_fun_jac_hess[i]);
    }
{%- endif %}


{% if dims.ns > 0 %}
    double Zl[NS];
    double Zu[NS];
    double zl[NS];
    double zu[NS];
    {% for j in range(end=dims.ns) %}
    Zl[{{ j }}] = {{ cost.Zl[j] }};
    {%- endfor %}

    {% for j in range(end=dims.ns) %}
    Zu[{{ j }}] = {{ cost.Zu[j] }};
    {%- endfor %}

    {% for j in range(end=dims.ns) %}
    zl[{{ j }}] = {{ cost.zl[j] }};
    {%- endfor %}

    {% for j in range(end=dims.ns) %}
    zu[{{ j }}] = {{ cost.zu[j] }};
    {%- endfor %}

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Zl", Zl);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Zu", Zu);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "zl", zl);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "zu", zu);
    }
{% endif %}

    // terminal cost
{% if cost.cost_type_e == "LINEAR_LS" or cost.cost_type_e == "NONLINEAR_LS" %}
{% if dims.ny_e > 0 %}
    double yref_e[NYN];
    {% for j in range(end=dims.ny_e) %}
    yref_e[{{ j }}] = {{ cost.yref_e[j] }};
    {%- endfor %}
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "yref", yref_e);

    double W_e[NYN*NYN];
    {% for j in range(end=dims.ny_e) %}
        {%- for k in range(end=dims.ny_e) %}
    W_e[{{ j }}+(NYN) * {{ k }}] = {{ cost.W_e[j][k] }};
        {%- endfor %}
    {%- endfor %}
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "W", W_e);

    {%- if cost.cost_type_e == "LINEAR_LS" %}
    double Vx_e[NYN*NX];
    {% for j in range(end=dims.ny_e) %}
        {%- for k in range(end=dims.nx) %}
    Vx_e[{{ j }}+(NYN) * {{ k }}] = {{ cost.Vx_e[j][k] }};
        {%- endfor %}
    {%- endfor %}
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "Vx", Vx_e);
    {%- endif %}

    {%- if cost.cost_type_e == "NONLINEAR_LS" %}
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "nls_y_fun", &cost_y_e_fun);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "nls_y_fun_jac", &cost_y_e_fun_jac_ut_xt);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "nls_y_hess", &cost_y_e_hess);
    {%- endif %}
{%- endif %}{# ny_e > 0 #}

{%- elif cost.cost_type_e == "EXTERNAL" %}
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "ext_cost_fun", &ext_cost_e_fun);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "ext_cost_fun_jac", &ext_cost_e_fun_jac);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "ext_cost_fun_jac_hess", &ext_cost_e_fun_jac_hess);
{%- endif %}

{% if dims.ns_e > 0 %}
    double Zl_e[NSN];
    double Zu_e[NSN];
    double zl_e[NSN];
    double zu_e[NSN];

    {% for j in range(end=dims.ns_e) %}
    Zl_e[{{ j }}] = {{ cost.Zl_e[j] }};
    {%- endfor %}

    {% for j in range(end=dims.ns_e) %}
    Zu_e[{{ j }}] = {{ cost.Zu_e[j] }};
    {%- endfor %}

    {% for j in range(end=dims.ns_e) %}
    zl_e[{{ j }}] = {{ cost.zl_e[j] }};
    {%- endfor %}

    {% for j in range(end=dims.ns_e) %}
    zu_e[{{ j }}] = {{ cost.zu_e[j] }};
    {%- endfor %}

    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "Zl", Zl_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "Zu", Zu_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "zl", zl_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "zu", zu_e);
{%- endif %}

    /**** Constraints ****/

    // bounds for initial stage
{% if dims.nbx_0 > 0 %}
    // x0
    int idxbx0[{{ dims.nbx_0 }}];
    {% for i in range(end=dims.nbx_0) %}
    idxbx0[{{ i }}] = {{ constraints.idxbx_0[i] }};
    {%- endfor %}

    double lbx0[{{ dims.nbx_0 }}];
    double ubx0[{{ dims.nbx_0 }}];
    {% for i in range(end=dims.nbx_0) %}
    lbx0[{{ i }}] = {{ constraints.lbx_0[i] }};
    ubx0[{{ i }}] = {{ constraints.ubx_0[i] }};
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxbx", idxbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", lbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", ubx0);
{% endif %}
{% if dims.nbxe_0 > 0 %}
    // idxbxe_0
    int idxbxe_0[{{ dims.nbxe_0 }}];
    {% for i in range(end=dims.nbxe_0) %}
    idxbxe_0[{{ i }}] = {{ constraints.idxbxe_0[i] }};
    {%- endfor %}
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxbxe", idxbxe_0);
{% endif %}

    /* constraints that are the same for initial and intermediate */
{%- if dims.nsbx > 0 %}
{# TODO: introduce nsbx0 & REMOVE SETTING lsbx, usbx for stage 0!!! move this block down!! #}
    // ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxsbx", idxsbx);
    // ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lsbx", lsbx);
    // ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "usbx", usbx);

    // soft bounds on x
    int idxsbx[NSBX];
    {% for i in range(end=dims.nsbx) %}
    idxsbx[{{ i }}] = {{ constraints.idxsbx[i] }};
    {%- endfor %}
    double lsbx[NSBX];
    double usbx[NSBX];
    {% for i in range(end=dims.nsbx) %}
    lsbx[{{ i }}] = {{ constraints.lsbx[i] }};
    usbx[{{ i }}] = {{ constraints.usbx[i] }};
    {%- endfor %}

    for (int i = 1; i < N; i++)
    {       
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxsbx", idxsbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lsbx", lsbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "usbx", usbx);
    }
{%- endif %}


{% if dims.nbu > 0 %}
    // u
    int idxbu[NBU];
    {% for i in range(end=dims.nbu) %}
    idxbu[{{ i }}] = {{ constraints.idxbu[i] }};
    {%- endfor %}
    double lbu[NBU];
    double ubu[NBU];
    {% for i in range(end=dims.nbu) %}
    lbu[{{ i }}] = {{ constraints.lbu[i] }};
    ubu[{{ i }}] = {{ constraints.ubu[i] }};
    {%- endfor %}

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxbu", idxbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lbu", lbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "ubu", ubu);
    }
{% endif %}

{% if dims.nsbu > 0 %}
    // set up soft bounds for u
    int idxsbu[NSBU];
    {% for i in range(end=dims.nsbu) %}
    idxsbu[{{ i }}] = {{ constraints.idxsbu[i] }};
    {%- endfor %}
    double lsbu[NSBU];
    double usbu[NSBU];
    {% for i in range(end=dims.nsbu) %}
    lsbu[{{ i }}] = {{ constraints.lsbu[i] }};
    usbu[{{ i }}] = {{ constraints.usbu[i] }};
    {%- endfor %}
    for (int i = 0; i < N; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxsbu", idxsbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lsbu", lsbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "usbu", usbu);
    }
{% endif %}

{% if dims.nsg > 0 %}
    // set up soft bounds for general linear constraints
    int idxsg[NSG];
    {% for i in range(end=dims.nsg) %}
    idxsg[{{ i }}] = {{ constraints.idxsg[i] }};
    {%- endfor %}
    double lsg[NSG];
    double usg[NSG];
    {% for i in range(end=dims.nsg) %}
    lsg[{{ i }}] = {{ constraints.lsg[i] }};
    usg[{{ i }}] = {{ constraints.usg[i] }};
    {%- endfor %}

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxsg", idxsg);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lsg", lsg);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "usg", usg);
    }
{% endif %}

{% if dims.nsh > 0 %}
    // set up soft bounds for nonlinear constraints
    int idxsh[NSH];
    {% for i in range(end=dims.nsh) %}
    idxsh[{{ i }}] = {{ constraints.idxsh[i] }};
    {%- endfor %}
    double lsh[NSH];
    double ush[NSH];
    {% for i in range(end=dims.nsh) %}
    lsh[{{ i }}] = {{ constraints.lsh[i] }};
    ush[{{ i }}] = {{ constraints.ush[i] }};
    {%- endfor %}

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxsh", idxsh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lsh", lsh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "ush", ush);
    }
{% endif %}

{% if dims.nsphi > 0 %}
    // set up soft bounds for convex-over-nonlinear constraints
    int idxsphi[NSPHI];
    {% for i in range(end=dims.nsphi) %}
    idxsphi[{{ i }}] = {{ constraints.idxsphi[i] }};
    {%- endfor %}
    double lsphi[NSPHI];
    double usphi[NSPHI];
    {% for i in range(end=dims.nsphi) %}
    lsphi[{{ i }}] = {{ constraints.lsphi[i] }};
    usphi[{{ i }}] = {{ constraints.usphi[i] }};
    {%- endfor %}

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxsphi", idxsphi);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lsphi", lsphi);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "usphi", usphi);
    }
{% endif %}

{% if dims.nbx > 0 %}
    // x
    int idxbx[NBX];
    {% for i in range(end=dims.nbx) %}
    idxbx[{{ i }}] = {{ constraints.idxbx[i] }};
    {%- endfor %}
    double lbx[NBX];
    double ubx[NBX];
    {% for i in range(end=dims.nbx) %}
    lbx[{{ i }}] = {{ constraints.lbx[i] }};
    ubx[{{ i }}] = {{ constraints.ubx[i] }};
    {%- endfor %}

    for (int i = 1; i < N; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxbx", idxbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lbx", lbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "ubx", ubx);
    }
{% endif %}

{% if dims.ng > 0 %}
    // set up general constraints for stage 0 to N-1 
    double D[NG*NU];
    double C[NG*NX];
    double lg[NG];
    double ug[NG];

    {% for j in range(end=dims.ng) %}
        {%- for k in range(end=dims.nu) %}
    D[{{ j }}+NG * {{ k }}] = {{ constraints.D[j][k] }};
        {%- endfor %}
    {%- endfor %}

    {% for j in range(end=dims.ng) %}
        {%- for k in range(end=dims.nx) %}
    C[{{ j }}+NG * {{ k }}] = {{ constraints.C[j][k] }};
        {%- endfor %}
    {%- endfor %}

    {% for i in range(end=dims.ng) %}
    lg[{{ i }}] = {{ constraints.lg[i] }};
    {%- endfor %}

    {% for i in range(end=dims.ng) %}
    ug[{{ i }}] = {{ constraints.ug[i] }};
    {%- endfor %}

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "D", D);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "C", C);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lg", lg);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "ug", ug);
    }
{% endif %}

{% if dims.nh > 0 %}
    // set up nonlinear constraints for stage 0 to N-1 
    double lh[NH];
    double uh[NH];

    {% for i in range(end=dims.nh) %}
    lh[{{ i }}] = {{ constraints.lh[i] }};
    {%- endfor %}

    {% for i in range(end=dims.nh) %}
    uh[{{ i }}] = {{ constraints.uh[i] }};
    {%- endfor %}
    
    for (int i = 0; i < N; i++)
    {
        // nonlinear constraints for stages 0 to N-1
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "nl_constr_h_fun_jac",
                                     &nl_constr_h_fun_jac[i]);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "nl_constr_h_fun",
                                    &nl_constr_h_fun[i]);
        {% if solver_options.hessian_approx == "EXACT" %}
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i,
                                      "nl_constr_h_fun_jac_hess", &nl_constr_h_fun_jac_hess[i]);
        {% endif %}
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lh", lh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "uh", uh);
    }
{% endif %}

{% if dims.nphi > 0 and constraints.constr_type == "BGP" %}
    // set up convex-over-nonlinear constraints for stage 0 to N-1 
    double lphi[NPHI];
    double uphi[NPHI];

    {% for i in range(end=dims.nphi) %}
    lphi[{{ i }}] = {{ constraints.lphi[i] }};
    {%- endfor %}

    {% for i in range(end=dims.nphi) %}
    uphi[{{ i }}] = {{ constraints.uphi[i] }};
    {%- endfor %}

    for (int i = 0; i < N; i++)
    {
        // ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "nl_constr_r_fun_jac", &r_constraint[i]);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i,
                      "nl_constr_phi_o_r_fun_phi_jac_ux_z_phi_hess_r_jac_ux", &phi_constraint[i]);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lphi", lphi);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "uphi", uphi);
    }
{% endif %}

    /* terminal constraints */
{% if dims.nbx_e > 0 %}
    // set up bounds for last stage
    // x
    int idxbx_e[NBXN];
    {% for i in range(end=dims.nbx_e) %}
    idxbx_e[{{ i }}] = {{ constraints.idxbx_e[i] }};
    {%- endfor %}
    double lbx_e[NBXN];
    double ubx_e[NBXN];
    {% for i in range(end=dims.nbx_e) %}
    lbx_e[{{ i }}] = {{ constraints.lbx_e[i] }};
    ubx_e[{{ i }}] = {{ constraints.ubx_e[i] }};
    {%- endfor %}
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "idxbx", idxbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lbx", lbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "ubx", ubx_e);
{%- endif %}

{% if dims.nsg_e > 0 %}
    // set up soft bounds for general linear constraints
    int idxsg_e[NSGN];
    {% for i in range(end=dims.nsg_e) %}
    idxsg_e[{{ i }}] = {{ constraints.idxsg_e[i] }};
    {%- endfor %}
    double lsg_e[NSGN];
    double usg_e[NSGN];
    {% for i in range(end=dims.nsg_e) %}
    lsg_e[{{ i }}] = {{ constraints.lsg_e[i] }};
    usg_e[{{ i }}] = {{ constraints.usg_e[i] }};
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "idxsg", idxsg_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lsg", lsg_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "usg", usg_e);
{%- endif %}

{% if dims.nsh_e > 0 %}
    // set up soft bounds for nonlinear constraints
    int idxsh_e[NSHN];
    {% for i in range(end=dims.nsh_e) %}
    idxsh_e[{{ i }}] = {{ constraints.idxsh_e[i] }};
    {%- endfor %}
    double lsh_e[NSHN];
    double ush_e[NSHN];
    {% for i in range(end=dims.nsh_e) %}
    lsh_e[{{ i }}] = {{ constraints.lsh_e[i] }};
    ush_e[{{ i }}] = {{ constraints.ush_e[i] }};
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "idxsh", idxsh_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lsh", lsh_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "ush", ush_e);
{%- endif %}

{% if dims.nsphi_e > 0 %}
    // set up soft bounds for convex-over-nonlinear constraints
    int idxsphi_e[NSPHIN];
    {% for i in range(end=dims.nsphi_e) %}
    idxsphi_e[{{ i }}] = {{ constraints.idxsphi_e[i] }};
    {%- endfor %}
    double lsphi_e[NSPHIN];
    double usphi_e[NSPHIN];
    {% for i in range(end=dims.nsphi_e) %}
    lsphi_e[{{ i }}] = {{ constraints.lsphi_e[i] }};
    usphi_e[{{ i }}] = {{ constraints.usphi_e[i] }};
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "idxsphi", idxsphi_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lsphi", lsphi_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "usphi", usphi_e);
{%- endif %}

{% if dims.nsbx_e > 0 %}
    // soft bounds on x
    int idxsbx_e[NSBXN];
    {% for i in range(end=dims.nsbx_e) %}
    idxsbx_e[{{ i }}] = {{ constraints.idxsbx_e[i] }};
    {%- endfor %}
    double lsbx_e[NSBXN];
    double usbx_e[NSBXN];
    {% for i in range(end=dims.nsbx_e) %}
    lsbx_e[{{ i }}] = {{ constraints.lsbx_e[i] }};
    usbx_e[{{ i }}] = {{ constraints.usbx_e[i] }};
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "idxsbx", idxsbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lsbx", lsbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "usbx", usbx_e);
{% endif %}

{% if dims.ng_e > 0 %}
    // set up general constraints for last stage 
    double C_e[NGN*NX];
    double lg_e[NGN];
    double ug_e[NGN];

    {% for j in range(end=dims.ng) %}
        {%- for k in range(end=dims.nx) %}
    C_e[{{ j }}+NG * {{ k }}] = {{ constraints.C_e[j][k] }};
        {%- endfor %}
    {%- endfor %}

    {% for i in range(end=dims.ng_e) %}
    lg_e[{{ i }}] = {{ constraints.lg_e[i] }};
    ug_e[{{ i }}] = {{ constraints.ug_e[i] }};
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "C", C_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lg", lg_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "ug", ug_e);
{%- endif %}

{% if dims.nh_e > 0 %}
    // set up nonlinear constraints for last stage 
    double lh_e[NHN];
    double uh_e[NHN];

    {% for i in range(end=dims.nh_e) %}
    lh_e[{{ i }}] = {{ constraints.lh_e[i] }};
    {%- endfor %}

    {% for i in range(end=dims.nh_e) %}
    uh_e[{{ i }}] = {{ constraints.uh_e[i] }};
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "nl_constr_h_fun_jac", &nl_constr_h_e_fun_jac);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "nl_constr_h_fun", &nl_constr_h_e_fun);
    {% if solver_options.hessian_approx == "EXACT" %}
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "nl_constr_h_fun_jac_hess", &nl_constr_h_e_fun_jac_hess);
    {% endif %}
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lh", lh_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "uh", uh_e);
{%- endif %}

{% if dims.nphi_e > 0 and constraints.constr_type_e == "BGP" %}
    // set up convex-over-nonlinear constraints for last stage 
    double lphi_e[NPHIN];
    double uphi_e[NPHIN];

    {% for i in range(end=dims.nphi_e) %}
    lphi_e[{{ i }}] = {{ constraints.lphi_e[i] }};
    {%- endfor %}

    {% for i in range(end=dims.nphi_e) %}
    uphi_e[{{ i }}] = {{ constraints.uphi_e[i] }};
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lphi", lphi_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "uphi", uphi_e);
    // ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "nl_constr_r_fun_jac", &r_e_constraint);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N,
                       "nl_constr_phi_o_r_fun_phi_jac_ux_z_phi_hess_r_jac_ux", &phi_e_constraint);
{% endif %}


    /************************************************
    *  opts
    ************************************************/

    nlp_opts = ocp_nlp_solver_opts_create(nlp_config, nlp_dims);

{% if solver_options.hessian_approx == "EXACT" %}
    bool nlp_solver_exact_hessian = true;
    // TODO: this if should not be needed! however, calling the setter with false leads to weird behavior. Investigate!
    if (nlp_solver_exact_hessian)
    {
        ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "exact_hess", &nlp_solver_exact_hessian);
    }
    int exact_hess_dyn = {{ solver_options.exact_hess_dyn }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "exact_hess_dyn", &exact_hess_dyn);

    int exact_hess_cost = {{ solver_options.exact_hess_cost }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "exact_hess_cost", &exact_hess_cost);

    int exact_hess_constr = {{ solver_options.exact_hess_constr }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "exact_hess_constr", &exact_hess_constr);
{%- endif -%}

{%- if dims.nz > 0 %}
    // TODO: these options are lower level -> should be encapsulated! maybe through hessian approx option.
    bool output_z_val = true;
    bool sens_algebraic_val = true;

    for (int i = 0; i < N; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_output_z", &output_z_val);
    for (int i = 0; i < N; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_sens_algebraic", &sens_algebraic_val);
{%- endif %}

    int num_steps_val = {{ solver_options.sim_method_num_steps }};
    for (int i = 0; i < N; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_num_steps", &num_steps_val);

    int ns_val = {{ solver_options.sim_method_num_stages }};
    for (int i = 0; i < N; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_num_stages", &ns_val);

    int newton_iter_val = {{ solver_options.sim_method_newton_iter }};
    for (int i = 0; i < N; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_newton_iter", &newton_iter_val);

    bool tmp_bool = {{ solver_options.sim_method_jac_reuse }};
    for (int i = 0; i < N; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_jac_reuse", &tmp_bool);

    double nlp_solver_step_length = {{ solver_options.nlp_solver_step_length }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "step_length", &nlp_solver_step_length);

    double levenberg_marquardt = {{ solver_options.levenberg_marquardt }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "levenberg_marquardt", &levenberg_marquardt);

    /* options QP solver */
{%- if solver_options.qp_solver is starting_with("PARTIAL_CONDENSING") %}
    int qp_solver_cond_N;

    {%- if solver_options.qp_solver_cond_N %}
    qp_solver_cond_N = {{ solver_options.qp_solver_cond_N }};
    {% else %}
    // NOTE: there is no condensing happening here!
    qp_solver_cond_N = N;
    {%- endif %}
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_cond_N", &qp_solver_cond_N);
{% endif %}

    int qp_solver_iter_max = {{ solver_options.qp_solver_iter_max }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_iter_max", &qp_solver_iter_max);

    {%- if solver_options.qp_solver_tol_stat %}
    double qp_solver_tol_stat = {{ solver_options.qp_solver_tol_stat }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_tol_stat", &qp_solver_tol_stat);
    {%- endif -%}

    {%- if solver_options.qp_solver_tol_eq %}
    double qp_solver_tol_eq = {{ solver_options.qp_solver_tol_eq }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_tol_eq", &qp_solver_tol_eq);
    {%- endif -%}

    {%- if solver_options.qp_solver_tol_ineq %}
    double qp_solver_tol_ineq = {{ solver_options.qp_solver_tol_ineq }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_tol_ineq", &qp_solver_tol_ineq);
    {%- endif -%}

    {%- if solver_options.qp_solver_tol_comp %}
    double qp_solver_tol_comp = {{ solver_options.qp_solver_tol_comp }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_tol_comp", &qp_solver_tol_comp);
    {%- endif -%}

{% if solver_options.nlp_solver_type == "SQP" %}
    // set SQP specific options
    double nlp_solver_tol_stat = {{ solver_options.nlp_solver_tol_stat }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_stat", &nlp_solver_tol_stat);

    double nlp_solver_tol_eq = {{ solver_options.nlp_solver_tol_eq }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_eq", &nlp_solver_tol_eq);

    double nlp_solver_tol_ineq = {{ solver_options.nlp_solver_tol_ineq }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_ineq", &nlp_solver_tol_ineq);

    double nlp_solver_tol_comp = {{ solver_options.nlp_solver_tol_comp }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_comp", &nlp_solver_tol_comp);

    int nlp_solver_max_iter = {{ solver_options.nlp_solver_max_iter }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "max_iter", &nlp_solver_max_iter);

    int initialize_t_slacks = {{ solver_options.initialize_t_slacks }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "initialize_t_slacks", &initialize_t_slacks);
{%- endif %}

    int print_level = {{ solver_options.print_level }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "print_level", &print_level);


    int ext_cost_num_hess = {{ solver_options.ext_cost_num_hess }};
{%- if cost.cost_type == "EXTERNAL" %}
    for (int i = 0; i < N; i++)
    {
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "cost_numerical_hessian", &ext_cost_num_hess);
    }
{%- endif %}
{%- if cost.cost_type_e == "EXTERNAL" %}
    ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, N, "cost_numerical_hessian", &ext_cost_num_hess);
{%- endif %}


    /* out */
    nlp_out = ocp_nlp_out_create(nlp_config, nlp_dims);

    // initialize primal solution
    double x0[{{ dims.nx }}];
{% if dims.nbx_0 == dims.nx %}
    // initialize with x0
    {% for item in constraints.lbx_0 %}
    x0[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}
{% else %}
    // initialize with zeros
    {% for i in range(end=dims.nx) %}
    x0[{{ i }}] = 0.0;
    {%- endfor %}
{%- endif %}

    double u0[NU];
    {% for i in range(end=dims.nu) %}
    u0[{{ i }}] = 0.0;
    {%- endfor %}

    for (int i = 0; i < N; i++)
    {
        // x0
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "x", x0);
        // u0
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "u", u0);
    }
    ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, N, "x", x0);
    
    nlp_solver = ocp_nlp_solver_create(nlp_config, nlp_dims, nlp_opts);


{% if dims.np > 0 %}
    // initialize parameters to nominal value
    double p[{{ dims.np }}];
    {% for i in range(end=dims.np) %}
    p[{{ i }}] = {{ parameter_values[i] }};
    {%- endfor %}

{% if solver_options.integrator_type == "IRK" %}
    for (int ii = 0; ii < N; ii++)
    {
        impl_dae_fun[ii].set_param(impl_dae_fun+ii, p);
        impl_dae_fun_jac_x_xdot_z[ii].set_param(impl_dae_fun_jac_x_xdot_z+ii, p);
        impl_dae_jac_x_xdot_u_z[ii].set_param(impl_dae_jac_x_xdot_u_z+ii, p);
        {%- if solver_options.hessian_approx == "EXACT" %}
        impl_dae_hess[ii].set_param(impl_dae_hess+ii, p);
        {%- endif %}
    }
{% elif solver_options.integrator_type == "ERK" %}
    for (int ii = 0; ii < N; ii++)
    {
        forw_vde_casadi[ii].set_param(forw_vde_casadi+ii, p);
        expl_ode_fun[ii].set_param(expl_ode_fun+ii, p);
    }
{% elif solver_options.integrator_type == "GNSF" %}
    for (int ii = 0; ii < N; ii++)
    {
        gnsf_phi_fun[ii].set_param(gnsf_phi_fun+ii, p);
        gnsf_phi_fun_jac_y[ii].set_param(gnsf_phi_fun_jac_y+ii, p);
        gnsf_phi_jac_y_uhat[ii].set_param(gnsf_phi_jac_y_uhat+ii, p);
        gnsf_f_lo_jac_x1_x1dot_u_z[ii].set_param(gnsf_f_lo_jac_x1_x1dot_u_z+ii, p);
    }
{% endif %}

    // cost
{%- if cost.cost_type == "NONLINEAR_LS" %}
    for (int ii = 0; ii < N; ii++)
    {
        cost_y_fun[ii].set_param(cost_y_fun+ii, p);
        cost_y_fun_jac_ut_xt[ii].set_param(cost_y_fun_jac_ut_xt+ii, p);
        cost_y_hess[ii].set_param(cost_y_hess+ii, p);
    }
{%- elif cost.cost_type == "EXTERNAL" %}
    for (int ii = 0; ii < N; ii++)
    {
        ext_cost_fun[ii].set_param(ext_cost_fun+ii, p);
        ext_cost_fun_jac[ii].set_param(ext_cost_fun_jac+ii, p);
        ext_cost_fun_jac_hess[ii].set_param(ext_cost_fun_jac_hess+ii, p);
    }
{%- endif %}

{%- if cost.cost_type_e == "NONLINEAR_LS" %}
    cost_y_e_fun.set_param(&cost_y_e_fun, p);
    cost_y_e_fun_jac_ut_xt.set_param(&cost_y_e_fun_jac_ut_xt, p);
    cost_y_e_hess.set_param(&cost_y_e_hess, p);
{%- elif cost.cost_type_e == "EXTERNAL" %}
    ext_cost_e_fun.set_param(&ext_cost_e_fun, p);
    ext_cost_e_fun_jac.set_param(&ext_cost_e_fun_jac, p);
    ext_cost_e_fun_jac_hess.set_param(&ext_cost_e_fun_jac_hess, p);
{%- endif %}

    // constraints
{%- if constraints.constr_type == "BGP" %}
    for (int ii = 0; ii < N; ii++)
    {
        // r_constraint[ii].set_param(r_constraint+ii, p);
        phi_constraint[ii].set_param(phi_constraint+ii, p);
    }
{%- elif dims.nh > 0 and constraints.constr_type == "BGH" %}

    for (int ii = 0; ii < N; ii++)
    {
        nl_constr_h_fun_jac[ii].set_param(nl_constr_h_fun_jac+ii, p);
        nl_constr_h_fun[ii].set_param(nl_constr_h_fun+ii, p);
    }
{%- if solver_options.hessian_approx == "EXACT" %}
    for (int ii = 0; ii < N; ii++)
    {
        nl_constr_h_fun_jac_hess[ii].set_param(nl_constr_h_fun_jac_hess+ii, p);
    }
{%- endif %}
{%- endif %}

{%- if constraints.constr_type_e == "BGP" %}
    // r_e_constraint.set_param(&r_e_constraint, p);
    phi_e_constraint.set_param(&phi_e_constraint, p);
{%- elif constraints.constr_type_e == "BGH" and dims.nh_e > 0 %}
    nl_constr_h_e_fun_jac.set_param(&nl_constr_h_e_fun_jac, p);
    nl_constr_h_e_fun.set_param(&nl_constr_h_e_fun, p);
{%- if solver_options.hessian_approx == "EXACT" %}
    nl_constr_h_e_fun_jac_hess.set_param(&nl_constr_h_e_fun_jac_hess, p);
{%- endif %}
{%- endif %}

{%- endif %}{# if dims.np #}

    status = ocp_nlp_precompute(nlp_solver, nlp_in, nlp_out);

    if (status != ACADOS_SUCCESS)
    {
        printf("\nocp_precompute failed!\n\n");
        exit(1);
    }

    return status;
}


int acados_update_params(int stage, double *p, int np)
{
    int solver_status = 0;

    int casadi_np = {{ dims.np }};
    if (casadi_np != np) {
        printf("acados_update_params: trying to set %i parameters for external functions."
            " External function has %i parameters. Exiting.\n", np, casadi_np);
        exit(1);
    }

{%- if dims.np > 0 %}
    if (stage < {{ dims.N }})
    {
    {%- if solver_options.integrator_type == "IRK" %}
        impl_dae_fun[stage].set_param(impl_dae_fun+stage, p);
        impl_dae_fun_jac_x_xdot_z[stage].set_param(impl_dae_fun_jac_x_xdot_z+stage, p);
        impl_dae_jac_x_xdot_u_z[stage].set_param(impl_dae_jac_x_xdot_u_z+stage, p);

        {%- if solver_options.hessian_approx == "EXACT" %}
        impl_dae_hess[stage].set_param(impl_dae_hess+stage, p);
        {%- endif %}
    {% elif solver_options.integrator_type == "ERK" %}
        forw_vde_casadi[stage].set_param(forw_vde_casadi+stage, p);
        expl_ode_fun[stage].set_param(expl_ode_fun+stage, p);

        {%- if solver_options.hessian_approx == "EXACT" %}
        hess_vde_casadi[stage].set_param(hess_vde_casadi+stage, p);
        {%- endif %}
    {% elif solver_options.integrator_type == "GNSF" %}
        gnsf_phi_fun[stage].set_param(gnsf_phi_fun+stage, p);
        gnsf_phi_fun_jac_y[stage].set_param(gnsf_phi_fun_jac_y+stage, p);
        gnsf_phi_jac_y_uhat[stage].set_param(gnsf_phi_jac_y_uhat+stage, p);

        gnsf_f_lo_jac_x1_x1dot_u_z[stage].set_param(gnsf_f_lo_jac_x1_x1dot_u_z+stage, p);
    {%- endif %}{# integrator_type #}

        // constraints
    {% if constraints.constr_type == "BGP" %}
        // r_constraint[stage].set_param(r_constraint+stage, p);
        phi_constraint[stage].set_param(phi_constraint+stage, p);
    {% elif constraints.constr_type == "BGH" and dims.nh > 0 %}
        nl_constr_h_fun_jac[stage].set_param(nl_constr_h_fun_jac+stage, p);
        nl_constr_h_fun[stage].set_param(nl_constr_h_fun+stage, p);
    {%- if solver_options.hessian_approx == "EXACT" %}
        nl_constr_h_fun_jac_hess[stage].set_param(nl_constr_h_fun_jac_hess+stage, p);
    {%- endif %}
    {%- endif %}

        // cost
    {%- if cost.cost_type == "NONLINEAR_LS" %}
        cost_y_fun[stage].set_param(cost_y_fun+stage, p);
        cost_y_fun_jac_ut_xt[stage].set_param(cost_y_fun_jac_ut_xt+stage, p);
        cost_y_hess[stage].set_param(cost_y_hess+stage, p);
    {%- elif cost.cost_type == "EXTERNAL" %}
        ext_cost_fun[stage].set_param(ext_cost_fun+stage, p);
        ext_cost_fun_jac[stage].set_param(ext_cost_fun_jac+stage, p);
        ext_cost_fun_jac_hess[stage].set_param(ext_cost_fun_jac_hess+stage, p);
    {%- endif %}

    }
    else // stage == N
    {
        // terminal shooting node has no dynamics
        // cost
    {%- if cost.cost_type_e == "NONLINEAR_LS" %}
        cost_y_e_fun.set_param(&cost_y_e_fun, p);
        cost_y_e_fun_jac_ut_xt.set_param(&cost_y_e_fun_jac_ut_xt, p);
        cost_y_e_hess.set_param(&cost_y_e_hess, p);
    {%- elif cost.cost_type_e == "EXTERNAL" %}
        ext_cost_e_fun.set_param(&ext_cost_e_fun, p);
        ext_cost_e_fun_jac_hess.set_param(&ext_cost_e_fun_jac_hess, p);
    {% endif %}
        // constraints
    {% if constraints.constr_type_e == "BGP" %}
        // r_e_constraint.set_param(&r_e_constraint, p);
        phi_e_constraint.set_param(&phi_e_constraint, p);
    {% elif constraints.constr_type_e == "BGH" and dims.nh_e > 0 %}
        nl_constr_h_e_fun_jac.set_param(&nl_constr_h_e_fun_jac, p);
        nl_constr_h_e_fun.set_param(&nl_constr_h_e_fun, p);
    {%- if solver_options.hessian_approx == "EXACT" %}
        nl_constr_h_e_fun_jac_hess[stage].set_param(nl_constr_h_e_fun_jac_hess+stage, p);
    {%- endif %}
    {% endif %}
    }
{% endif %}{# if dims.np #}

    return solver_status;
}



int acados_solve()
{
    // solve NLP 
    int solver_status = ocp_nlp_solve(nlp_solver, nlp_in, nlp_out);

    return solver_status;
}


int acados_free()
{
    // free memory
    ocp_nlp_solver_opts_destroy(nlp_opts);
    ocp_nlp_in_destroy(nlp_in);
    ocp_nlp_out_destroy(nlp_out);
    ocp_nlp_solver_destroy(nlp_solver);
    ocp_nlp_dims_destroy(nlp_dims);
    ocp_nlp_config_destroy(nlp_config);
    ocp_nlp_plan_destroy(nlp_solver_plan);

    /* free external function */
    // dynamics
{%- if solver_options.integrator_type == "IRK" %}
    for (int i = 0; i < {{ dims.N }}; i++)
    {
        external_function_param_casadi_free(&impl_dae_fun[i]);
        external_function_param_casadi_free(&impl_dae_fun_jac_x_xdot_z[i]);
        external_function_param_casadi_free(&impl_dae_jac_x_xdot_u_z[i]);
    {%- if solver_options.hessian_approx == "EXACT" %}
        external_function_param_casadi_free(&impl_dae_hess[i]);
    {%- endif %}
    }
    free(impl_dae_fun);
    free(impl_dae_fun_jac_x_xdot_z);
    free(impl_dae_jac_x_xdot_u_z);
    {%- if solver_options.hessian_approx == "EXACT" %}
    free(impl_dae_hess);
    {%- endif %}

{%- elif solver_options.integrator_type == "ERK" %}
    for (int i = 0; i < {{ dims.N }}; i++)
    {
        external_function_param_casadi_free(&forw_vde_casadi[i]);
        external_function_param_casadi_free(&expl_ode_fun[i]);
    {%- if solver_options.hessian_approx == "EXACT" %}
        external_function_param_casadi_free(&hess_vde_casadi[i]);
    {%- endif %}
    }
    free(forw_vde_casadi);
    free(expl_ode_fun);
    {%- if solver_options.hessian_approx == "EXACT" %}
    free(hess_vde_casadi);
    {%- endif %}

{%- elif solver_options.integrator_type == "GNSF" %}
    for (int i = 0; i < {{ dims.N }}; i++)
    {
        external_function_param_casadi_free(&gnsf_phi_fun[i]);
        external_function_param_casadi_free(&gnsf_phi_fun_jac_y[i]);
        external_function_param_casadi_free(&gnsf_phi_jac_y_uhat[i]);
        external_function_param_casadi_free(&gnsf_f_lo_jac_x1_x1dot_u_z[i]);
        external_function_param_casadi_free(&gnsf_get_matrices_fun[i]);
    }
    free(gnsf_phi_fun);
    free(gnsf_phi_fun_jac_y);
    free(gnsf_phi_jac_y_uhat);
    free(gnsf_f_lo_jac_x1_x1dot_u_z);
    free(gnsf_get_matrices_fun);
{%- endif %}

    // cost
{%- if cost.cost_type == "NONLINEAR_LS" %}
    for (int i = 0; i < {{ dims.N }}; i++)
    {
        external_function_param_casadi_free(&cost_y_fun[i]);
        external_function_param_casadi_free(&cost_y_fun_jac_ut_xt[i]);
        external_function_param_casadi_free(&cost_y_hess[i]);
    }
    free(cost_y_fun);
    free(cost_y_fun_jac_ut_xt);
    free(cost_y_hess);
{%- elif cost.cost_type == "EXTERNAL" %}
    for (int i = 0; i < {{ dims.N }}; i++)
    {
        external_function_param_casadi_free(&ext_cost_fun[i]);
        external_function_param_casadi_free(&ext_cost_fun_jac[i]);
        external_function_param_casadi_free(&ext_cost_fun_jac_hess[i]);
    }
    free(ext_cost_fun);
    free(ext_cost_fun_jac_hess);
{%- endif %}
{%- if cost.cost_type_e == "NONLINEAR_LS" %}
    external_function_param_casadi_free(&cost_y_e_fun);
    external_function_param_casadi_free(&cost_y_e_fun_jac_ut_xt);
    external_function_param_casadi_free(&cost_y_e_hess);
{%- elif cost.cost_type_e == "EXTERNAL" %}
    external_function_param_casadi_free(&ext_cost_e_fun);
    external_function_param_casadi_free(&ext_cost_e_fun_jac_hess);
{%- endif %}

    // constraints
{%- if constraints.constr_type == "BGH" and dims.nh > 0 %}
    for (int i = 0; i < {{ dims.N }}; i++)
    {
        external_function_param_casadi_free(&nl_constr_h_fun_jac[i]);
        external_function_param_casadi_free(&nl_constr_h_fun[i]);
    }
  {%- if solver_options.hessian_approx == "EXACT" %}
    for (int i = 0; i < {{ dims.N }}; i++)
    {
        external_function_param_casadi_free(&nl_constr_h_fun_jac_hess[i]);
    }
  {%- endif %}
    free(nl_constr_h_fun_jac);
    free(nl_constr_h_fun);
  {%- if solver_options.hessian_approx == "EXACT" %}
    free(nl_constr_h_fun_jac_hess);
  {%- endif %}

{%- elif constraints.constr_type == "BGP" and dims.nphi > 0 %}
    for (int i = 0; i < {{ dims.N }}; i++)
    {
        external_function_param_casadi_free(&phi_constraint[i]);
    }
    free(phi_constraint);
{%- endif %}

{%- if constraints.constr_type_e == "BGH" and dims.nh_e > 0 %}
    external_function_param_casadi_free(&nl_constr_h_e_fun_jac);
    external_function_param_casadi_free(&nl_constr_h_e_fun);
{%- if solver_options.hessian_approx == "EXACT" %}
    external_function_param_casadi_free(&nl_constr_h_e_fun_jac_hess);
{%- endif %}
{%- elif constraints.constr_type_e == "BGP" and dims.nphi_e > 0 %}
    external_function_param_casadi_free(&phi_e_constraint);
{%- endif %}

    return 0;
}

ocp_nlp_in * acados_get_nlp_in() { return  nlp_in; }
ocp_nlp_out * acados_get_nlp_out() { return  nlp_out; }
ocp_nlp_solver * acados_get_nlp_solver() { return  nlp_solver; }
ocp_nlp_config * acados_get_nlp_config() { return  nlp_config; }
void * acados_get_nlp_opts() { return  nlp_opts; }
ocp_nlp_dims * acados_get_nlp_dims() { return  nlp_dims; }
ocp_nlp_plan * acados_get_nlp_plan() { return  nlp_solver_plan; }


void acados_print_stats()
{
    int sqp_iter, stat_m, stat_n, tmp_int;
    ocp_nlp_get(nlp_config, nlp_solver, "sqp_iter", &sqp_iter);
    ocp_nlp_get(nlp_config, nlp_solver, "stat_n", &stat_n);
    ocp_nlp_get(nlp_config, nlp_solver, "stat_m", &stat_m);

    {% set stat_n_max = 10 %}
    double stat[{{ solver_options.nlp_solver_max_iter * stat_n_max }}];
    ocp_nlp_get(nlp_config, nlp_solver, "statistics", stat);

    int nrow = sqp_iter+1 < stat_m ? sqp_iter+1 : stat_m;

    printf("iter\tres_stat\tres_eq\t\tres_ineq\tres_comp\tqp_stat\tqp_iter\n");
    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < stat_n + 1; j++)
        {
            if (j == 0 || j > 4)
            {
                tmp_int = (int) stat[i + j * nrow];
                printf("%d\t", tmp_int);
            }
            else
            {
                printf("%e\t", stat[i + j * nrow]);
            }
        }
        printf("\n");
    }
}
