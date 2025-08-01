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

{% set dims_0 = phases_dims | first %}
{% set cost_0 = cost | first %}
{% set constraints_0 = constraints | first %}
{% set model_0 = model | first %}


{% set cost_e = cost | last %}
{% set constraints_e = constraints | last %}
{% set dims_e = phases_dims | last %}
{% set model_e = model | last %}


{%- set nx_values = [] -%}
{%- for jj in range(end=n_phases) %}
    {%- set_global nx_values = nx_values | concat(with=(phases_dims[jj].nx)) %}
{%- endfor %}
{%- set nx_max = nx_values | sort | last %}

{%- set nu_values = [] -%}
{%- for jj in range(end=n_phases) %}
    {%- set_global nu_values = nu_values | concat(with=(phases_dims[jj].nu)) %}
{%- endfor %}
{%- set nu_max = nu_values | sort | last %}

{%- set np_values = [] -%}
{%- for jj in range(end=n_phases) %}
    {%- set_global np_values = np_values | concat(with=(phases_dims[jj].np)) %}
{%- endfor %}
{%- set np_max = np_values | sort | last %}


// standard
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
// acados
// #include "acados/utils/print.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

// example specific
    {%- for jj in range(end=n_phases) %}{# phases loop !#}
#include "{{ model[jj].name }}_model/{{ model[jj].name }}_model.h"

{%- if phases_dims[jj].nh > 0 or phases_dims[jj].nh_e > 0 or phases_dims[jj].nh_0 > 0 or phases_dims[jj].nphi > 0 or phases_dims[jj].nphi_e > 0 or phases_dims[jj].nphi_0 > 0 %}
#include "{{ model[jj].name }}_constraints/{{ model[jj].name }}_constraints.h"
{%- endif %}
{%- if cost[jj].cost_type != "LINEAR_LS" or cost[jj].cost_type_e != "LINEAR_LS" or cost[jj].cost_type_0 != "LINEAR_LS" %}
#include "{{ model[jj].name }}_cost/{{ model[jj].name }}_cost.h"
{%- endif %}
    {%- endfor %}{# for jj in range(end=n_phases) #}

{% if phases_dims[0].n_global_data > 0 %}
#include "{{ name }}_p_global_precompute_fun.h"
{%- endif %}

{%- if not solver_options.custom_update_filename %}
    {%- set custom_update_filename = "" %}
{% else %}
    {%- set custom_update_filename = solver_options.custom_update_filename %}
{%- endif %}
{%- if not solver_options.custom_update_header_filename %}
    {%- set custom_update_header_filename = "" %}
{% else %}
    {%- set custom_update_header_filename = solver_options.custom_update_header_filename %}
{%- endif %}
{%- if custom_update_header_filename != "" %}
#include "{{ custom_update_header_filename }}"
{%- endif %}

#include "acados_solver_{{ name }}.h"


#define {{ name | upper }}_N      {{ N_horizon }}

{%- for jj in range(end=n_phases) %}
#define NP_{{ jj }}     {{ phases_dims[jj].np }}
{%- endfor %}


{{ name }}_solver_capsule * {{ name }}_acados_create_capsule(void)
{
    void* capsule_mem = malloc(sizeof({{ name }}_solver_capsule));
    {{ name }}_solver_capsule *capsule = ({{ name }}_solver_capsule *) capsule_mem;

    return capsule;
}


int {{ name }}_acados_free_capsule({{ name }}_solver_capsule *capsule)
{
    free(capsule);
    return 0;
}


int {{ name }}_acados_create({{ name }}_solver_capsule* capsule)
{
    int N_shooting_intervals = {{ name | upper }}_N;
    double* new_time_steps = NULL; // NULL -> don't alter the code generated time-steps
    return {{ name }}_acados_create_with_discretization(capsule, N_shooting_intervals, new_time_steps);
}


/**
 * Internal function for {{ name }}_acados_create: step 1
 */
void {{ name }}_acados_create_set_plan(ocp_nlp_plan_t* nlp_solver_plan, const int N)
{
    assert(N == nlp_solver_plan->N);

    /************************************************
    *  plan
    ************************************************/

    nlp_solver_plan->nlp_solver = {{ solver_options.nlp_solver_type }};

    nlp_solver_plan->ocp_qp_solver_plan.qp_solver = {{ solver_options.qp_solver }};
    nlp_solver_plan->relaxed_ocp_qp_solver_plan.qp_solver = {{ solver_options.qp_solver }};

    nlp_solver_plan->regularization = {{ solver_options.regularize_method }};
    nlp_solver_plan->globalization = {{ solver_options.globalization }};

    nlp_solver_plan->nlp_cost[0] = {{ cost[0].cost_type_0 }};
    nlp_solver_plan->nlp_constraints[0] = {{ constraints[0].constr_type_0 }};

{%- for jj in range(end=n_phases) %}{# phases loop !#}
    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        nlp_solver_plan->nlp_cost[i] = {{ cost[jj].cost_type }};
        nlp_solver_plan->nlp_constraints[i] = {{ constraints[jj].constr_type }};
    }
    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
      {%- if mocp_opts.integrator_type[jj] == "DISCRETE" %}
        nlp_solver_plan->nlp_dynamics[i] = DISCRETE_MODEL;
        // discrete dynamics does not need sim solver option, this field is ignored
        nlp_solver_plan->sim_solver_plan[i].sim_solver = INVALID_SIM_SOLVER;
      {%- else %}
        nlp_solver_plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
        nlp_solver_plan->sim_solver_plan[i].sim_solver = {{ mocp_opts.integrator_type[jj] }};
      {%- endif %}
    }
{%- endfor %}

    nlp_solver_plan->nlp_cost[N] = {{ cost_e.cost_type_e }};
    nlp_solver_plan->nlp_constraints[N] = {{ constraints_e.constr_type_e }};
}



/**
 * Internal function for {{ name }}_acados_create: step 2
 */
ocp_nlp_dims* {{ name }}_acados_create_setup_dimensions({{ name }}_solver_capsule* capsule)
{
    ocp_nlp_plan_t* nlp_solver_plan = capsule->nlp_solver_plan;
    const int N = nlp_solver_plan->N;
    ocp_nlp_config* nlp_config = capsule->nlp_config;
    int i;

    /************************************************
    *  dimensions
    ************************************************/
    #define NINTNP1MEMS 18
    int* intNp1mem = (int*)malloc( (N+1)*sizeof(int)*NINTNP1MEMS );

    int* nx    = intNp1mem + (N+1)*0;
    int* nu    = intNp1mem + (N+1)*1;
    int* nbx   = intNp1mem + (N+1)*2;
    int* nbu   = intNp1mem + (N+1)*3;
    int* nsbx  = intNp1mem + (N+1)*4;
    int* nsbu  = intNp1mem + (N+1)*5;
    int* nsg   = intNp1mem + (N+1)*6;
    int* nsh   = intNp1mem + (N+1)*7;
    int* nsphi = intNp1mem + (N+1)*8;
    int* ns    = intNp1mem + (N+1)*9;
    int* ng    = intNp1mem + (N+1)*10;
    int* nh    = intNp1mem + (N+1)*11;
    int* nphi  = intNp1mem + (N+1)*12;
    int* nz    = intNp1mem + (N+1)*13;
    int* ny    = intNp1mem + (N+1)*14;
    int* nr    = intNp1mem + (N+1)*15;
    int* nbxe  = intNp1mem + (N+1)*16;
    int* np    = intNp1mem + (N+1)*17;

{%- for jj in range(end=n_phases) %}{# phases loop !#}
    for (i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        // common
        nx[i] = {{ phases_dims[jj].nx }};
        nu[i] = {{ phases_dims[jj].nu }};
        nz[i] = {{ phases_dims[jj].nz }};
        ns[i] = {{ phases_dims[jj].ns }};
        np[i] = {{ phases_dims[jj].np }};
        // cost
        ny[i] = {{ phases_dims[jj].ny }};
        // constraints
        nbu[i] = {{ phases_dims[jj].nbu }};
        nbx[i] = {{ phases_dims[jj].nbx }};
        ng[i] = {{ phases_dims[jj].ng }};
        nh[i] = {{ phases_dims[jj].nh }};
        nphi[i] = {{ phases_dims[jj].nphi }};
        nr[i] = {{ phases_dims[jj].nr }};
        // slacks
        nsbu[i] = {{ phases_dims[jj].nsbu }};
        nsbx[i] = {{ phases_dims[jj].nsbx }};
        nsg[i] = {{ phases_dims[jj].nsg }};
        nsh[i] = {{ phases_dims[jj].nsh }};
        nsphi[i] = {{ phases_dims[jj].nsphi }};
        nbxe[i] = 0;
    }
{%- endfor %}{# phases loop !#}

    /* initial node*/
    i = 0;
    // common
    nx[i] = {{ phases_dims[0].nx }};
    nu[i] = {{ phases_dims[0].nu }};
    nz[i] = {{ phases_dims[0].nz }};
    ns[i] = {{ phases_dims[0].ns_0 }};
    np[i] = {{ phases_dims[0].np }};
    // cost
    ny[i] = {{ phases_dims[0].ny_0 }};
    // constraints
    nbu[i] = {{ phases_dims[0].nbu }};
    nbx[i] = {{ phases_dims[0].nbx_0 }};
    nbxe[i] = {{ phases_dims[0].nbxe_0 }};
    ng[i] = {{ phases_dims[0].ng }};
    nh[i] = {{ phases_dims[0].nh_0 }};
    nphi[i] = {{ phases_dims[0].nphi_0 }};
    nr[i] = {{ phases_dims[0].nr_0 }};
    // slacks
    nsbu[i] = {{ phases_dims[0].nsbu }};
    nsbx[i] = {{ 0 }};
    nsg[i] = {{ phases_dims[0].nsg }};
    nsh[i] = {{ phases_dims[0].nsh_0 }};
    nsphi[i] = {{ phases_dims[0].nsphi_0 }};

    /* terminal node */
    // common
    i = N;
    nx[i] = {{ dims_e.nx }};
    nu[i] = {{ 0 }};
    nz[i] = {{ dims_e.nz }};
    ns[i] = {{ dims_e.ns_e }};
    np[i] = {{ dims_e.np }};
    // cost
    ny[i] = {{ dims_e.ny_e }};
    // constraints
    nbu[i] = {{ 0 }};
    nbx[i] = {{ dims_e.nbx_e }};
    ng[i] = {{ dims_e.ng_e }};
    nh[i] = {{ dims_e.nh_e }};
    nphi[i] = {{ dims_e.nphi_e }};
    nr[i] = {{ dims_e.nr_e }};
    // slacks
    nsbu[i] = {{ 0 }};
    nsbx[i] = {{ dims_e.nsbx_e }};
    nsg[i] = {{ dims_e.nsg_e }};
    nsh[i] = {{ dims_e.nsh_e }};
    nsphi[i] = {{ dims_e.nsphi_e }};
    nbxe[i] = 0;

    /* create and set ocp_nlp_dims */
    ocp_nlp_dims * nlp_dims = ocp_nlp_dims_create(nlp_config);
{# options that are always set #}
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nx", nx);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nu", nu);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nz", nz);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "ns", ns);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "np", np);

    ocp_nlp_dims_set_global(nlp_config, nlp_dims, "np_global", {{ dims_0.np_global }});
    ocp_nlp_dims_set_global(nlp_config, nlp_dims, "n_global_data", {{ dims_0.n_global_data }});

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

{# set module depedent dimension -- INITIAL NODE #}
{%- if cost[0].cost_type_0 == "NONLINEAR_LS" or cost[0].cost_type_0 == "LINEAR_LS" or cost[0].cost_type_0 == "CONVEX_OVER_NONLINEAR"%}
    ocp_nlp_dims_set_cost(nlp_config, nlp_dims, 0, "ny", &ny[0]);
{%- endif %}

{%- if constraints[0].constr_type_0 == "BGH" %}
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, 0, "nh", &nh[0]);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, 0, "nsh", &nsh[0]);
{%- elif constraints[0].constr_type_0 == "BGP" %}
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, 0, "nr", &nr[0]);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, 0, "nphi", &nphi[0]);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, 0, "nsphi", &nsphi[0]);
{%- endif %}

{# set module depedent dimension -- PATH #}
{%- for jj in range(end=n_phases) %}{# phases loop !#}
    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
      {%- if cost[jj].cost_type == "NONLINEAR_LS" or cost[jj].cost_type == "LINEAR_LS" or cost[jj].cost_type == "CONVEX_OVER_NONLINEAR"%}
        ocp_nlp_dims_set_cost(nlp_config, nlp_dims, i, "ny", &ny[i]);
      {%- endif %}
      {%- if constraints[jj].constr_type == "BGH" %}
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nh", &nh[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsh", &nsh[i]);
      {%- elif constraints[jj].constr_type == "BGP" %}
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nr", &nr[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nphi", &nphi[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsphi", &nsphi[i]);
      {%- endif %}
    }

  {%- if mocp_opts.integrator_type[jj] == "GNSF" -%}
    // GNSF specific dimensions
    int gnsf_nx1 = {{ phases_dims[jj].gnsf_nx1 }};
    int gnsf_nz1 = {{ phases_dims[jj].gnsf_nz1 }};
    int gnsf_nout = {{ phases_dims[jj].gnsf_nout }};
    int gnsf_ny = {{ phases_dims[jj].gnsf_ny }};
    int gnsf_nuhat = {{ phases_dims[jj].gnsf_nuhat }};

    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_dims_set_dynamics(nlp_config, nlp_dims, i, "gnsf_nx1", &gnsf_nx1);
        ocp_nlp_dims_set_dynamics(nlp_config, nlp_dims, i, "gnsf_nz1", &gnsf_nz1);
        ocp_nlp_dims_set_dynamics(nlp_config, nlp_dims, i, "gnsf_nout", &gnsf_nout);
        ocp_nlp_dims_set_dynamics(nlp_config, nlp_dims, i, "gnsf_ny", &gnsf_ny);
        ocp_nlp_dims_set_dynamics(nlp_config, nlp_dims, i, "gnsf_nuhat", &gnsf_nuhat);
    }
  {%- endif %}

  {%- if mocp_opts.cost_discretization[jj] == "INTEGRATOR" %}
    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
        ocp_nlp_dims_set_dynamics(nlp_config, nlp_dims, i, "ny", &ny[i]);
  {%- endif %}

{%- endfor %}{# phases loop !#}


{# set module depedent dimension -- TERMINAL #}
{%- if constraints_e.constr_type_e == "BGH" %}
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "nh", &nh[N]);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "nsh", &nsh[N]);
{%- elif constraints_e.constr_type_e == "BGP" %}
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "nr", &nr[N]);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "nphi", &nphi[N]);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "nsphi", &nsphi[N]);
{%- endif %}
{%- if cost_e.cost_type_e == "NONLINEAR_LS" or cost_e.cost_type_e == "LINEAR_LS" or cost_e.cost_type_e == "CONVEX_OVER_NONLINEAR"%}
    ocp_nlp_dims_set_cost(nlp_config, nlp_dims, N, "ny", &ny[N]);
{%- endif %}

    free(intNp1mem);

    return nlp_dims;
}



/**
 * Internal function for {{ name }}_acados_create: step 3
 */
void {{ name }}_acados_create_setup_functions({{ name }}_solver_capsule* capsule)
{
    const int N = capsule->nlp_solver_plan->N;

    /************************************************
    *  external functions
    ************************************************/

#define MAP_CASADI_FNC(__CAPSULE_FNC__, __MODEL_BASE_FNC__) do{ \
        capsule->__CAPSULE_FNC__.casadi_fun = & __MODEL_BASE_FNC__ ;\
        capsule->__CAPSULE_FNC__.casadi_n_in = & __MODEL_BASE_FNC__ ## _n_in; \
        capsule->__CAPSULE_FNC__.casadi_n_out = & __MODEL_BASE_FNC__ ## _n_out; \
        capsule->__CAPSULE_FNC__.casadi_sparsity_in = & __MODEL_BASE_FNC__ ## _sparsity_in; \
        capsule->__CAPSULE_FNC__.casadi_sparsity_out = & __MODEL_BASE_FNC__ ## _sparsity_out; \
        capsule->__CAPSULE_FNC__.casadi_work = & __MODEL_BASE_FNC__ ## _work; \
        external_function_external_param_casadi_create(&capsule->__CAPSULE_FNC__, &ext_fun_opts); \
    } while(false)

    external_function_opts ext_fun_opts;
    external_function_opts_set_to_default(&ext_fun_opts);

{% if phases_dims[0].n_global_data > 0 %}
    // NOTE: p_global_precompute_fun cannot use external_workspace!!!
    ext_fun_opts.external_workspace = false;
    capsule->p_global_precompute_fun.casadi_fun = &{{ name }}_p_global_precompute_fun;
    capsule->p_global_precompute_fun.casadi_work = &{{ name }}_p_global_precompute_fun_work;
    capsule->p_global_precompute_fun.casadi_sparsity_in = &{{ name }}_p_global_precompute_fun_sparsity_in;
    capsule->p_global_precompute_fun.casadi_sparsity_out = &{{ name }}_p_global_precompute_fun_sparsity_out;
    capsule->p_global_precompute_fun.casadi_n_in = &{{ name }}_p_global_precompute_fun_n_in;
    capsule->p_global_precompute_fun.casadi_n_out = &{{ name }}_p_global_precompute_fun_n_out;
    external_function_casadi_create(&capsule->p_global_precompute_fun, &ext_fun_opts);

    ext_fun_opts.with_global_data = true;
{%- endif %}

    ext_fun_opts.external_workspace = true;

{# INITIAL #}
{%- if constraints[0].constr_type_0 == "BGH" and phases_dims[0].nh_0 > 0 %}
    MAP_CASADI_FNC(nl_constr_h_0_fun_jac, {{ model[0].name }}_constr_h_0_fun_jac_uxt_zt);
    MAP_CASADI_FNC(nl_constr_h_0_fun, {{ model[0].name }}_constr_h_0_fun);

    {%- if solver_options.hessian_approx == "EXACT" %}
    MAP_CASADI_FNC(nl_constr_h_0_fun_jac_hess, {{ model[0].name }}_constr_h_0_fun_jac_uxt_zt_hess);
    {% endif %}
    {%- if solver_options.with_solution_sens_wrt_params %}
    MAP_CASADI_FNC(nl_constr_h_0_jac_p_hess_xu_p, {{ model[0].name }}_constr_h_0_jac_p_hess_xu_p);
    {%- endif %}
    {%- if solver_options.with_value_sens_wrt_params %}
    MAP_CASADI_FNC(nl_constr_h_0_adj_p, {{ model[0].name }}_constr_h_0_adj_p);
    {%- endif %}
{%- elif constraints[0].constr_type_0 == "BGP" %}
    // convex-over-nonlinear constraint
    MAP_CASADI_FNC(phi_0_constraint_fun, {{ model[0].name }}_phi_0_constraint_fun);
    MAP_CASADI_FNC(phi_0_constraint_fun_jac_hess, {{ model[0].name }}_phi_0_constraint_fun_jac_hess);
{%- endif %}


{%- if cost[0].cost_type_0 == "NONLINEAR_LS" %}
    // nonlinear least squares function
    MAP_CASADI_FNC(cost_y_0_fun, {{ model[0].name }}_cost_y_0_fun);
    MAP_CASADI_FNC(cost_y_0_fun_jac_ut_xt, {{ model[0].name }}_cost_y_0_fun_jac_ut_xt);
    {%- if solver_options.hessian_approx == "EXACT" %}
    MAP_CASADI_FNC(cost_y_0_hess, {{ model[0].name }}_cost_y_0_hess);
    {%- endif %}

{%- elif cost[0].cost_type_0 == "CONVEX_OVER_NONLINEAR" %}
    // convex-over-nonlinear cost
    MAP_CASADI_FNC(conl_cost_0_fun, {{ model[0].name }}_conl_cost_0_fun);
    MAP_CASADI_FNC(conl_cost_0_fun_jac_hess, {{ model[0].name }}_conl_cost_0_fun_jac_hess);

{%- elif cost[0].cost_type_0 == "EXTERNAL" %}
    // external cost
    {%- if cost[0].cost_ext_fun_type_0 == "casadi" %}
    MAP_CASADI_FNC(ext_cost_0_fun, {{ model[0].name }}_cost_ext_cost_0_fun);
    {%- else %}
    capsule->ext_cost_0_fun.fun = &{{ cost[0].cost_function_ext_cost_0 }};
    external_function_external_param_{{ cost[0].cost_ext_fun_type_0 }}_create(&capsule->ext_cost_0_fun, &ext_fun_opts);
    {%- endif %}

    // external cost
    {%- if cost[0].cost_ext_fun_type_0 == "casadi" %}
    MAP_CASADI_FNC(ext_cost_0_fun_jac, {{ model[0].name }}_cost_ext_cost_0_fun_jac);
    {%- else %}
    capsule->ext_cost_0_fun_jac.fun = &{{ cost[0].cost_function_ext_cost_0 }};
    external_function_external_param_{{ cost[0].cost_ext_fun_type_0 }}_create(&capsule->ext_cost_0_fun_jac, &ext_fun_opts);
    {%- endif %}

    // external cost
    {%- if cost[0].cost_ext_fun_type_0 == "casadi" %}
    MAP_CASADI_FNC(ext_cost_0_fun_jac_hess, {{ model[0].name }}_cost_ext_cost_0_fun_jac_hess);
    {%- else %}
    capsule->ext_cost_0_fun_jac_hess.fun = &{{ cost[0].cost_function_ext_cost_0 }};
    external_function_external_param_{{ cost[0].cost_ext_fun_type_0 }}_create(&capsule->ext_cost_0_fun_jac_hess, &ext_fun_opts);
    {%- endif %}
{%- endif %}



/////////////// PATH
    int n_path, n_cost_path;

    {%- for jj in range(end=n_phases) %}{# phases loop !#}
    n_path = {{ end_idx[jj] - start_idx[jj]}};
    n_cost_path = {{ end_idx[jj] - cost_start_idx[jj] }};
{%- if constraints[jj].constr_type == "BGH" and phases_dims[jj].nh > 0  %}
    capsule->nl_constr_h_fun_jac_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++) {
        MAP_CASADI_FNC(nl_constr_h_fun_jac_{{ jj }}[i], {{ model[jj].name }}_constr_h_fun_jac_uxt_zt);
    }
    capsule->nl_constr_h_fun_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++) {
        MAP_CASADI_FNC(nl_constr_h_fun_{{ jj }}[i], {{ model[jj].name }}_constr_h_fun);
    }
    {% if solver_options.hessian_approx == "EXACT" %}
    capsule->nl_constr_h_fun_jac_hess_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++) {
        MAP_CASADI_FNC(nl_constr_h_fun_jac_hess_{{ jj }}[i], {{ model[jj].name }}_constr_h_fun_jac_uxt_zt_hess);
    }
    {% endif %}
    {%- if solver_options.with_solution_sens_wrt_params %}
    capsule->nl_constr_h_jac_p_hess_xu_p_{ jj } = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*(N-1));
    for (int i = 0; i < n_cost_path; i++) {
        MAP_CASADI_FNC(nl_constr_h_jac_p_hess_xu_p_{{ jj }}[i], {{ model[jj].name }}_constr_h_jac_p_hess_xu_p);
    }
    {%- endif %}
    {%- if solver_options.with_value_sens_wrt_params %}
    capsule->nl_constr_h_adj_p_{ jj } = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*(N-1));
    for (int i = 0; i < n_cost_path; i++) {
        MAP_CASADI_FNC(nl_constr_h_adj_p_{{ jj }}[i], {{ model[jj].name }}_constr_h_adj_p);
    }
    {%- endif %}
{% elif constraints[jj].constr_type == "BGP" %}
    capsule->phi_constraint_fun_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    capsule->phi_constraint_fun_jac_hess_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        // convex-over-nonlinear constraint
        MAP_CASADI_FNC(phi_constraint_fun_{{ jj }}[i], {{ model[jj].name }}_phi_constraint_fun);
        MAP_CASADI_FNC(phi_constraint_fun_jac_hess_{{ jj }}[i], {{ model[jj].name }}_phi_constraint_fun_jac_hess);
    }
{%- endif %}

{% if mocp_opts.integrator_type[jj] == "ERK" %}
    // explicit ode
    capsule->expl_vde_forw_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(expl_vde_forw_{{ jj }}[i], {{ model[jj].name }}_expl_vde_forw);
    }

    capsule->expl_vde_adj_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(expl_vde_adj_{{ jj }}[i], {{ model[jj].name }}_expl_vde_adj);
    }

    capsule->expl_ode_fun_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(expl_ode_fun_{{ jj }}[i], {{ model[jj].name }}_expl_ode_fun);
    }

    {%- if solver_options.hessian_approx == "EXACT" %}
    capsule->expl_ode_hess_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(expl_ode_hess_{{ jj }}[i], {{ model[jj].name }}_expl_ode_hess);
    }
    {%- endif %}

{% elif mocp_opts.integrator_type[jj] == "IRK" %}
    // implicit dae
    capsule->impl_dae_fun_{{ jj }} = (external_function_external_param_{{ model[jj].dyn_ext_fun_type }} *) malloc(sizeof(external_function_external_param_{{ model[jj].dyn_ext_fun_type }})*n_path);
    for (int i = 0; i < n_path; i++) {
    {%- if model[jj].dyn_ext_fun_type == "casadi" %}
        MAP_CASADI_FNC(impl_dae_fun_{{ jj }}[i], {{ model[jj].name }}_impl_dae_fun);
    {%- else %}
        capsule->impl_dae_fun_{{ jj }}[i].fun = &{{ model[jj].dyn_impl_dae_fun }};
        external_function_external_param_{{ model[jj].dyn_ext_fun_type }}_create(&capsule->impl_dae_fun_{{ jj }}[i], &ext_fun_opts);
    {%- endif %}
    }

    capsule->impl_dae_fun_jac_x_xdot_z_{{ jj }} = (external_function_external_param_{{ model[jj].dyn_ext_fun_type }} *) malloc(sizeof(external_function_external_param_{{ model[jj].dyn_ext_fun_type }})*n_path);
    for (int i = 0; i < n_path; i++) {
    {%- if model[jj].dyn_ext_fun_type == "casadi" %}
        MAP_CASADI_FNC(impl_dae_fun_jac_x_xdot_z_{{ jj }}[i], {{ model[jj].name }}_impl_dae_fun_jac_x_xdot_z);
    {%- else %}
        capsule->impl_dae_fun_jac_x_xdot_z_{{ jj }}[i].fun = &{{ model[jj].dyn_impl_dae_fun_jac }};
        external_function_external_param_{{ model[jj].dyn_ext_fun_type }}_create(&capsule->impl_dae_fun_jac_x_xdot_z_{{ jj }}[i], &ext_fun_opts);
    {%- endif %}
    }

    capsule->impl_dae_jac_x_xdot_u_z_{{ jj }} = (external_function_external_param_{{ model[jj].dyn_ext_fun_type }} *) malloc(sizeof(external_function_external_param_{{ model[jj].dyn_ext_fun_type }})*n_path);
    for (int i = 0; i < n_path; i++) {
    {%- if model[jj].dyn_ext_fun_type == "casadi" %}
        MAP_CASADI_FNC(impl_dae_jac_x_xdot_u_z_{{ jj }}[i], {{ model[jj].name }}_impl_dae_jac_x_xdot_u_z);
    {%- else %}
        capsule->impl_dae_jac_x_xdot_u_z_{{ jj }}[i].fun = &{{ model[jj].dyn_impl_dae_jac }};
        external_function_external_param_{{ model[jj].dyn_ext_fun_type }}_create(&capsule->impl_dae_jac_x_xdot_u_z_{{ jj }}[i], &ext_fun_opts);
    {%- endif %}
    }

    {%- if solver_options.hessian_approx == "EXACT" %}
    capsule->impl_dae_hess_{{ jj }} = (external_function_external_param_{{ model[jj].dyn_ext_fun_type }} *) malloc(sizeof(external_function_external_param_{{ model[jj].dyn_ext_fun_type }})*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(impl_dae_hess_{{ jj }}[i], {{ model[jj].name }}_impl_dae_hess);
    }
    {%- endif %}
{% elif mocp_opts.integrator_type[jj] == "LIFTED_IRK" %}
    // external functions (implicit model)
    capsule->impl_dae_fun_{{ jj }} = (external_function_external_param_{{ model[jj].dyn_ext_fun_type }} *) malloc(sizeof(external_function_external_param_{{ model[jj].dyn_ext_fun_type }})*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(impl_dae_fun_{{ jj }}[i], {{ model[jj].name }}_impl_dae_fun);
    }

    capsule->impl_dae_fun_jac_x_xdot_u_{{ jj }} = (external_function_external_param_{{ model[jj].dyn_ext_fun_type }} *) malloc(sizeof(external_function_external_param_{{ model[jj].dyn_ext_fun_type }})*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(impl_dae_fun_jac_x_xdot_u_{{ jj }}[i], {{ model[jj].name }}_impl_dae_fun_jac_x_xdot_u);
    }

{% elif mocp_opts.integrator_type[jj] == "GNSF" %}
    {% if model[jj].gnsf_purely_linear != 1 %}
    capsule->gnsf_phi_fun_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(gnsf_phi_fun_{{ jj }}[i], {{ model[jj].name }}_gnsf_phi_fun);
    }

    capsule->gnsf_phi_fun_jac_y_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(gnsf_phi_fun_jac_y_{{ jj }}[i], {{ model[jj].name }}_gnsf_phi_fun_jac_y);
    }

    capsule->gnsf_phi_jac_y_uhat_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(gnsf_phi_jac_y_uhat_{{ jj }}[i], {{ model[jj].name }}_gnsf_phi_jac_y_uhat);
    }

    {% if model[jj].gnsf_nontrivial_f_LO == 1 %}
    capsule->gnsf_f_lo_jac_x1_x1dot_u_z_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(gnsf_f_lo_jac_x1_x1dot_u_z_{{ jj }}[i], {{ model[jj].name }}_gnsf_f_lo_fun_jac_x1k1uz);
    }
    {%- endif %}
    {%- endif %}
    capsule->gnsf_get_matrices_fun_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(gnsf_get_matrices_fun_{{ jj }}[i], {{ model[jj].name }}_gnsf_get_matrices_fun);
    }
{% elif mocp_opts.integrator_type[jj] == "DISCRETE" %}
    // discrete dynamics
    capsule->discr_dyn_phi_fun_{{ jj }} = (external_function_external_param_{{ model[jj].dyn_ext_fun_type }} *) malloc(sizeof(external_function_external_param_{{ model[jj].dyn_ext_fun_type }})*n_path);
    for (int i = 0; i < n_path; i++)
    {
        {%- if model[jj].dyn_ext_fun_type == "casadi" %}
        MAP_CASADI_FNC(discr_dyn_phi_fun_{{ jj }}[i], {{ model[jj].name }}_dyn_disc_phi_fun);
        {%- else %}
        capsule->discr_dyn_phi_fun_{{ jj }}[i].fun = &{{ model[jj].dyn_disc_fun }};
        external_function_external_param_{{ model[jj].dyn_ext_fun_type }}_create(&capsule->discr_dyn_phi_fun_{{ jj }}[i], &ext_fun_opts);
        {%- endif %}
    }

    capsule->discr_dyn_phi_fun_jac_ut_xt_{{ jj }} = (external_function_external_param_{{ model[jj].dyn_ext_fun_type }} *) malloc(sizeof(external_function_external_param_{{ model[jj].dyn_ext_fun_type }})*n_path);
    for (int i = 0; i < n_path; i++)
    {
        {%- if model[jj].dyn_ext_fun_type == "casadi" %}
        MAP_CASADI_FNC(discr_dyn_phi_fun_jac_ut_xt_{{ jj }}[i], {{ model[jj].name }}_dyn_disc_phi_fun_jac);
        {%- else %}
        capsule->discr_dyn_phi_fun_jac_ut_xt_{{ jj }}[i].fun = &{{ model[jj].dyn_disc_fun_jac }};
        external_function_external_param_{{ model[jj].dyn_ext_fun_type }}_create(&capsule->discr_dyn_phi_fun_jac_ut_xt_{{ jj }}[i], &ext_fun_opts);
        {%- endif %}
    }

  {%- if solver_options.hessian_approx == "EXACT" %}
    capsule->discr_dyn_phi_fun_jac_ut_xt_hess_{{ jj }} = (external_function_external_param_{{ model[jj].dyn_ext_fun_type }} *) malloc(sizeof(external_function_external_param_{{ model[jj].dyn_ext_fun_type }})*n_path);
    for (int i = 0; i < n_path; i++)
    {
        {%- if model[jj].dyn_ext_fun_type == "casadi" %}
        MAP_CASADI_FNC(discr_dyn_phi_fun_jac_ut_xt_hess_{{ jj }}[i], {{ model[jj].name }}_dyn_disc_phi_fun_jac_hess);
        {%- else %}
        capsule->discr_dyn_phi_fun_jac_ut_xt_hess_{{ jj }}[i].fun = &{{ model[jj].dyn_disc_fun_jac_hess }};
        external_function_external_param_{{ model[jj].dyn_ext_fun_type }}_create(&capsule->discr_dyn_phi_fun_jac_ut_xt_hess_{{ jj }}[i], &ext_fun_opts);
        {%- endif %}
    }
  {%- endif %}
{%- endif %}



{%- if cost[jj].cost_type == "NONLINEAR_LS" %}
    // nonlinear least squares cost
    capsule->cost_y_fun_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        MAP_CASADI_FNC(cost_y_fun_{{ jj }}[i], {{ model[jj].name }}_cost_y_fun);
    }

    capsule->cost_y_fun_jac_ut_xt_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        MAP_CASADI_FNC(cost_y_fun_jac_ut_xt_{{ jj }}[i], {{ model[jj].name }}_cost_y_fun_jac_ut_xt);
    }

    {%- if solver_options.hessian_approx == "EXACT" %}
    capsule->cost_y_hess_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        MAP_CASADI_FNC(cost_y_hess_{{ jj }}[i], {{ model[jj].name }}_cost_y_hess);
    }
    {%- endif %}


{%- elif cost[jj].cost_type == "CONVEX_OVER_NONLINEAR" %}
    // convex-over-nonlinear cost
    capsule->conl_cost_fun_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        MAP_CASADI_FNC(conl_cost_fun_{{ jj }}[i], {{ model[jj].name }}_conl_cost_fun);
    }
    capsule->conl_cost_fun_jac_hess_{{ jj }} = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        MAP_CASADI_FNC(conl_cost_fun_jac_hess_{{ jj }}[i], {{ model[jj].name }}_conl_cost_fun_jac_hess);
    }

{%- elif cost[jj].cost_type == "EXTERNAL" %}
    // external cost
    capsule->ext_cost_fun_{{ jj }} = (external_function_external_param_{{ cost[jj].cost_ext_fun_type }} *) malloc(sizeof(external_function_external_param_{{ cost[jj].cost_ext_fun_type }})*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        {%- if cost[jj].cost_ext_fun_type == "casadi" %}
        MAP_CASADI_FNC(ext_cost_fun_{{ jj }}[i], {{ model[jj].name }}_cost_ext_cost_fun);
        {%- else %}
        capsule->ext_cost_fun_{{ jj }}[i].fun = &{{ cost[jj].cost_function_ext_cost }};
        external_function_external_param_{{ cost[jj].cost_ext_fun_type }}_create(&capsule->ext_cost_fun_{{ jj }}[i], &ext_fun_opts);
        {%- endif %}
    }

    capsule->ext_cost_fun_jac_{{ jj }} = (external_function_external_param_{{ cost[jj].cost_ext_fun_type }} *) malloc(sizeof(external_function_external_param_{{ cost[jj].cost_ext_fun_type }})*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        {%- if cost[jj].cost_ext_fun_type == "casadi" %}
        MAP_CASADI_FNC(ext_cost_fun_jac_{{ jj }}[i], {{ model[jj].name }}_cost_ext_cost_fun_jac);
        {%- else %}
        capsule->ext_cost_fun_jac_{{ jj }}[i].fun = &{{ cost[jj].cost_function_ext_cost }};
        external_function_external_param_{{ cost[jj].cost_ext_fun_type }}_create(&capsule->ext_cost_fun_jac_{{ jj }}[i], &ext_fun_opts);
        {%- endif %}
    }

    capsule->ext_cost_fun_jac_hess_{{ jj }} = (external_function_external_param_{{ cost[jj].cost_ext_fun_type }} *) malloc(sizeof(external_function_external_param_{{ cost[jj].cost_ext_fun_type }})*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        {%- if cost[jj].cost_ext_fun_type == "casadi" %}
        MAP_CASADI_FNC(ext_cost_fun_jac_hess_{{ jj }}[i], {{ model[jj].name }}_cost_ext_cost_fun_jac_hess);
        {%- else %}
        capsule->ext_cost_fun_jac_hess_{{ jj }}[i].fun = &{{ cost[jj].cost_function_ext_cost }};
        external_function_external_param_{{ cost[jj].cost_ext_fun_type }}_create(&capsule->ext_cost_fun_jac_hess_{{ jj }}[i], &ext_fun_opts);
        {%- endif %}
    }
{%- endif %}
{%- endfor %}



{# TERMINAL NODE #}
{%- if constraints_e.constr_type_e == "BGH" and dims_e.nh_e > 0 %}
    MAP_CASADI_FNC(nl_constr_h_e_fun_jac, {{ model_e.name }}_constr_h_e_fun_jac_uxt_zt);
    MAP_CASADI_FNC(nl_constr_h_e_fun, {{ model_e.name }}_constr_h_e_fun);

    {%- if solver_options.hessian_approx == "EXACT" %}
    MAP_CASADI_FNC(nl_constr_h_e_fun_jac_hess, {{ model_e.name }}_constr_h_e_fun_jac_uxt_zt_hess);
    {% endif %}
    {% if solver_options.with_solution_sens_wrt_params %}
    MAP_CASADI_FNC(nl_constr_h_e_jac_p_hess_xu_p, {{ model_e.name }}_constr_h_e_jac_p_hess_xu_p);
    {%- endif %}
    {% if solver_options.with_value_sens_wrt_params %}
    MAP_CASADI_FNC(nl_constr_h_e_adj_p, {{ model_e.name }}_constr_h_e_adj_p);
    {%- endif %}
{%- elif constraints_e.constr_type_e == "BGP" %}
    // convex-over-nonlinear constraint
    MAP_CASADI_FNC(phi_e_constraint_fun, {{ model_e.name }}_phi_e_constraint_fun);
    MAP_CASADI_FNC(phi_e_constraint_fun_jac_hess, {{ model_e.name }}_phi_e_constraint_fun_jac_hess);
{%- endif %}

{%- if cost_e.cost_type_e == "NONLINEAR_LS" %}
    // nonlinear least square function
    MAP_CASADI_FNC(cost_y_e_fun, {{ model_e.name }}_cost_y_e_fun);
    MAP_CASADI_FNC(cost_y_e_fun_jac_ut_xt, {{ model_e.name }}_cost_y_e_fun_jac_ut_xt);
    {%- if solver_options.hessian_approx == "EXACT" %}
    MAP_CASADI_FNC(cost_y_e_hess, {{ model_e.name }}_cost_y_e_hess);
    {%- endif %}
{%- elif cost_e.cost_type_e == "CONVEX_OVER_NONLINEAR" %}
    // convex-over-nonlinear cost
    MAP_CASADI_FNC(conl_cost_e_fun, {{ model_e.name }}_conl_cost_e_fun);
    MAP_CASADI_FNC(conl_cost_e_fun_jac_hess, {{ model_e.name }}_conl_cost_e_fun_jac_hess);

{%- elif cost_e.cost_type_e == "EXTERNAL" %}
    // external cost - function
    {%- if cost_e.cost_ext_fun_type_e == "casadi" %}
    MAP_CASADI_FNC(ext_cost_e_fun, {{ model_e.name }}_cost_ext_cost_e_fun);
    {%- else %}
    capsule->ext_cost_e_fun.fun = &{{ cost_e.cost_function_ext_cost_e }};
    external_function_external_param_{{ cost_e.cost_ext_fun_type_e }}_create(&capsule->ext_cost_e_fun, &ext_fun_opts);
    {%- endif %}

    // external cost - jacobian
    {%- if cost_e.cost_ext_fun_type_e == "casadi" %}
    MAP_CASADI_FNC(ext_cost_e_fun_jac, {{ model_e.name }}_cost_ext_cost_e_fun_jac);
    {%- else %}
    capsule->ext_cost_e_fun_jac.fun = &{{ cost_e.cost_function_ext_cost_e }};
    external_function_external_param_{{ cost_e.cost_ext_fun_type_e }}_create(&capsule->ext_cost_e_fun_jac, &ext_fun_opts);
    {%- endif %}

    // external cost - hessian
    {%- if cost_e.cost_ext_fun_type_e == "casadi" %}
    MAP_CASADI_FNC(ext_cost_e_fun_jac_hess, {{ model_e.name }}_cost_ext_cost_e_fun_jac_hess);
    {%- else %}
    capsule->ext_cost_e_fun_jac_hess.fun = &{{ cost_e.cost_function_ext_cost_e }};
    external_function_external_param_{{ cost_e.cost_ext_fun_type_e }}_create(&capsule->ext_cost_e_fun_jac_hess, &ext_fun_opts);
    {%- endif %}
{%- endif %}

#undef MAP_CASADI_FNC
}



/**
 * Internal function for {{ name }}_acados_create: step 4
 */
void {{ name }}_acados_create_set_default_parameters({{ name }}_solver_capsule* capsule) {

    double* p = calloc({{ np_max }}, sizeof(double));
    // initialize parameters to nominal value
{%- for jj in range(end=n_phases) %}{# phases loop !#}
    {%- for item in parameter_values[jj] %}
        {%- if item != 0 %}
    p[{{ loop.index0 }}] = {{ item }};
        {%- endif %}
    {%- endfor %}

    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++) {
        {{ name }}_acados_update_params(capsule, i, p, NP_{{ jj }});
    }
{%- endfor %}
    free(p);

{%- if phases_dims[0].np_global > 0 %}
    double* p_global = calloc({{ phases_dims[0].np_global }}, sizeof(double));
    // initialize global parameters to nominal value
    {%- for item in p_global_values %}
        {%- if item != 0 %}
    p_global[{{ loop.index0 }}] = {{ item }};
        {%- endif %}
    {%- endfor %}

    {{ name }}_acados_set_p_global_and_precompute_dependencies(capsule, p_global, {{ phases_dims[0].np_global }});

    free(p_global);
{%- endif %}
}




/**
 * Internal function for {{ name }}_acados_create: step 5
 */
void {{ name }}_acados_create_setup_nlp_in({{ name }}_solver_capsule* capsule, int N)
{
    assert(N == capsule->nlp_solver_plan->N);
    ocp_nlp_config* nlp_config = capsule->nlp_config;
    ocp_nlp_dims* nlp_dims = capsule->nlp_dims;

    int tmp_int = 0;

    /************************************************
    *  nlp_in
    ************************************************/
    ocp_nlp_in * nlp_in = capsule->nlp_in;
    /************************************************
    *  nlp_out
    ************************************************/
    ocp_nlp_out * nlp_out = capsule->nlp_out;


    // set up time_steps
    {%- set time_steps_all_equal = true -%}
    {%- set val = solver_options.time_steps[0] %}
    {%- for j in range(start=1, end=N_horizon) %}
        {%- if val != solver_options.time_steps[j] %}
            {%- set_global time_steps_all_equal = false %}
            {%- break %}
        {%- endif %}
    {%- endfor %}

{% if time_steps_all_equal == true %}{# all time_steps are identical #}
    double time_step = {{ solver_options.time_steps[0] }};
    for (int i = 0; i < N; i++)
    {
        ocp_nlp_in_set(nlp_config, nlp_dims, nlp_in, i, "Ts", &time_step);
    }
{%- else -%}{# time_steps are varying #}
    double* time_steps = malloc(N*sizeof(double));
    {%- for j in range(end=solver_options.N_horizon) %}
    time_steps[{{ j }}] = {{ solver_options.time_steps[j] }};
    {%- endfor %}
    for (int i = 0; i < N; i++)
    {
        ocp_nlp_in_set(nlp_config, nlp_dims, nlp_in, i, "Ts", &time_steps[i]);
    }
    free(time_steps);
{%- endif %}
    // set cost scaling
    double* cost_scaling = malloc((N+1)*sizeof(double));
    {%- for j in range(end=solver_options.N_horizon+1) %}
    cost_scaling[{{ j }}] = {{ solver_options.cost_scaling[j] }};
    {%- endfor %}
    for (int i = 0; i <= N; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "scaling", &cost_scaling[i]);
    }
    free(cost_scaling);

    /* INITIAL NODE */
{%- if dims_0.ny_0 != 0 %}
    double* yref_0 = calloc({{ dims_0.ny_0 }}, sizeof(double));
    // change only the non-zero elements:
    {%- for j in range(end=dims_0.ny_0) %}
        {%- if cost[0].yref_0[j] != 0 %}
    yref_0[{{ j }}] = {{ cost[0].yref_0[j] }};
        {%- endif %}
    {%- endfor %}
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, 0, "yref", yref_0);
    free(yref_0);
  {%- if cost[0].cost_type_0 == "NONLINEAR_LS" or cost[0].cost_type_0 == "LINEAR_LS" %}

    double* W_0 = calloc({{ dims_0.ny_0 }}*{{ dims_0.ny_0 }}, sizeof(double));
    // change only the non-zero elements:
    {%- for j in range(end=dims_0.ny_0) %}
        {%- for k in range(end=dims_0.ny_0) %}
            {%- if cost[0].W_0[j][k] != 0 %}
    W_0[{{ j }}+({{ dims_0.ny_0 }}) * {{ k }}] = {{ cost[0].W_0[j][k] }};
            {%- endif %}
        {%- endfor %}
    {%- endfor %}
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, 0, "W", W_0);
    free(W_0);
  {%- endif %}

  {%- if cost[0].cost_type_0 == "LINEAR_LS" %}
    double* Vx_0 = calloc({{ dims_0.ny_0 }}*{{ dims_0.nx }}, sizeof(double));
    // change only the non-zero elements:
    {%- for j in range(end=dims_0.ny_0) %}
        {%- for k in range(end=dims_0.nx) %}
            {%- if cost[0].Vx_0[j][k] != 0 %}
    Vx_0[{{ j }}+({{ dims_0.ny_0 }}) * {{ k }}] = {{ cost[0].Vx_0[j][k] }};
            {%- endif %}
        {%- endfor %}
    {%- endfor %}
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, 0, "Vx", Vx_0);
    free(Vx_0);

    {%- if dims_0.nu > 0 %}
    double* Vu_0 = calloc({{ dims_0.ny_0 }}*{{ dims_0.nu }}, sizeof(double));
    // change only the non-zero elements:
    {%- for j in range(end=dims_0.ny_0) %}
        {%- for k in range(end=dims_0.nu) %}
            {%- if cost[0].Vu_0[j][k] != 0 %}
    Vu_0[{{ j }}+({{ dims_0.ny_0 }}) * {{ k }}] = {{ cost[0].Vu_0[j][k] }};
            {%- endif %}
        {%- endfor %}
    {%- endfor %}
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, 0, "Vu", Vu_0);
    free(Vu_0);
    {%- endif %}

    {%- if dims_0.nz > 0 %}
    double* Vz_0 = calloc({{ dims_0.ny_0 }}*{{ dims_0.nz }}, sizeof(double));
    // change only the non-zero elements:
    {%- for j in range(end=dims_0.ny_0) %}
        {%- for k in range(end=dims_0.nz) %}
            {%- if cost[0].Vz_0[j][k] != 0 %}
    Vz_0[{{ j }}+({{ dims_0.ny_0 }}) * {{ k }}] = {{ cost[0].Vz_0[j][k] }};
            {%- endif %}
        {%- endfor %}
    {%- endfor %}
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, 0, "Vz", Vz_0);
    free(Vz_0);
    {%- endif %}{# nz > 0 #}
  {%- endif %}{# LINEAR_LS #}
{%- endif %}{# ny_0 != 0 #}


{%- if cost[0].cost_type_0 == "NONLINEAR_LS" %}
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "nls_y_fun", &capsule->cost_y_0_fun);
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "nls_y_fun_jac", &capsule->cost_y_0_fun_jac_ut_xt);
    {%- if solver_options.hessian_approx == "EXACT" %}
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "nls_y_hess", &capsule->cost_y_0_hess);
    {%- endif %}
{%- elif cost[0].cost_type_0 == "CONVEX_OVER_NONLINEAR" %}
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "conl_cost_fun", &capsule->conl_cost_0_fun);
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "conl_cost_fun_jac_hess", &capsule->conl_cost_0_fun_jac_hess);
{%- elif cost[0].cost_type_0 == "EXTERNAL" %}
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "ext_cost_fun", &capsule->ext_cost_0_fun);
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "ext_cost_fun_jac", &capsule->ext_cost_0_fun_jac);
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "ext_cost_fun_jac_hess", &capsule->ext_cost_0_fun_jac_hess);
{%- endif %}

{%- if mocp_opts.cost_discretization[0] == "INTEGRATOR" %}
  {%- if cost[0].cost_type_0 == "NONLINEAR_LS" %}
    ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "nls_y_fun_jac", &capsule->cost_y_0_fun_jac_ut_xt);
    ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "nls_y_fun", &capsule->cost_y_0_fun);
  {%- elif cost[0].cost_type_0 == "CONVEX_OVER_NONLINEAR" %}
    ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "conl_cost_fun", &capsule->conl_cost_0_fun);
    ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "conl_cost_fun_jac_hess", &capsule->conl_cost_0_fun_jac_hess);
  {%- endif %}
{%- endif %}


{% if dims_0.ns_0 > 0 %}
    // slacks terminal
    double* zlu0_mem = calloc(4*{{ dims_0.ns_0 }}, sizeof(double));
    double* Zl_0 = zlu0_mem+{{ dims_0.ns_0 }}*0;
    double* Zu_0 = zlu0_mem+{{ dims_0.ns_0 }}*1;
    double* zl_0 = zlu0_mem+{{ dims_0.ns_0 }}*2;
    double* zu_0 = zlu0_mem+{{ dims_0.ns_0 }}*3;

    // change only the non-zero elements:
    {%- for j in range(end=dims_0.ns_0) %}
        {%- if cost_0.Zl_0[j] != 0 %}
    Zl_0[{{ j }}] = {{ cost_0.Zl_0[j] }};
        {%- endif %}
    {%- endfor %}

    {%- for j in range(end=dims_0.ns_0) %}
        {%- if cost_0.Zu_0[j] != 0 %}
    Zu_0[{{ j }}] = {{ cost_0.Zu_0[j] }};
        {%- endif %}
    {%- endfor %}

    {%- for j in range(end=dims_0.ns_0) %}
        {%- if cost_0.zl_0[j] != 0 %}
    zl_0[{{ j }}] = {{ cost_0.zl_0[j] }};
        {%- endif %}
    {%- endfor %}

    {%- for j in range(end=dims_0.ns_0) %}
        {%- if cost_0.zu_0[j] != 0 %}
    zu_0[{{ j }}] = {{ cost_0.zu_0[j] }};
        {%- endif %}
    {%- endfor %}

    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, 0, "Zl", Zl_0);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, 0, "Zu", Zu_0);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, 0, "zl", zl_0);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, 0, "zu", zu_0);
    free(zlu0_mem);
{%- endif %}

    // constraints at initial node
{%- if dims_0.nbx_0 > 0 %}
    // x0
    int nbx_0 = {{ dims_0.nbx_0 }};
    int* idxbx0 = malloc(nbx_0 * sizeof(int));
    {%- for i in range(end=dims_0.nbx_0) %}
    idxbx0[{{ i }}] = {{ constraints_0.idxbx_0[i] }};
    {%- endfor %}

    double* lubx0 = calloc(2*nbx_0, sizeof(double));
    double* lbx0 = lubx0;
    double* ubx0 = lubx0 + nbx_0;
    // change only the non-zero elements:
    {%- for i in range(end=dims_0.nbx_0) %}
        {%- if constraints_0.lbx_0[i] != 0 %}
    lbx0[{{ i }}] = {{ constraints_0.lbx_0[i] }};
        {%- endif %}
        {%- if constraints_0.ubx_0[i] != 0 %}
    ubx0[{{ i }}] = {{ constraints_0.ubx_0[i] }};
        {%- endif %}
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "idxbx", idxbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "lbx", lbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "ubx", ubx0);
    free(idxbx0);
    free(lubx0);
{%- endif %}

{%- if dims_0.nbxe_0 > 0 %}
    // idxbxe_0
    int* idxbxe_0 = malloc({{ dims_0.nbxe_0 }} * sizeof(int));
    {%- for i in range(end=dims_0.nbxe_0) %}
    idxbxe_0[{{ i }}] = {{ constraints_0.idxbxe_0[i] }};
    {%- endfor %}
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "idxbxe", idxbxe_0);
    free(idxbxe_0);
{%- endif %}


{% if dims_0.nh_0 > 0 %}
    // set up nonlinear constraints for first stage
    double* luh_0 = calloc(2*{{ dims_0.nh_0 }}, sizeof(double));
    double* lh_0 = luh_0;
    double* uh_0 = luh_0 + {{ dims_0.nh_0 }};
    {%- for i in range(end=dims_0.nh_0) %}
        {%- if constraints_0.lh_0[i] != 0 %}
    lh_0[{{ i }}] = {{ constraints_0.lh_0[i] }};
        {%- endif %}
    {%- endfor %}

    {%- for i in range(end=dims_0.nh_0) %}
        {%- if constraints_0.uh_0[i] != 0 %}
    uh_0[{{ i }}] = {{ constraints_0.uh_0[i] }};
        {%- endif %}
    {%- endfor %}

    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "nl_constr_h_fun_jac", &capsule->nl_constr_h_0_fun_jac);
    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "nl_constr_h_fun", &capsule->nl_constr_h_0_fun);
    {% if solver_options.hessian_approx == "EXACT" %}
    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "nl_constr_h_fun_jac_hess",
                                  &capsule->nl_constr_h_0_fun_jac_hess);
    {% endif %}
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "lh", lh_0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "uh", uh_0);
    {% if solver_options.with_solution_sens_wrt_params %}
    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "nl_constr_h_jac_p_hess_xu_p",
                                  &capsule->nl_constr_h_0_jac_p_hess_xu_p);
    {% endif %}
    {% if solver_options.with_value_sens_wrt_params %}
    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "nl_constr_h_adj_p",
                                  &capsule->nl_constr_h_0_adj_p);
    {% endif %}
    free(luh_0);
{%- elif dims_0.nphi_0 > 0 and constraints_0.constr_type_0 == "BGP" %}
    // set up convex-over-nonlinear constraints for first stage
    double* luphi_0 = calloc(2*{{ dims_0.nphi_0 }}, sizeof(double));
    double* lphi_0 = luphi_0;
    double* uphi_0 = luphi_0 + {{ dims_0.nphi_0 }};
    {%- for i in range(end=dims_0.nphi_0) %}
        {%- if constraints_0.lphi_0[i] != 0 %}
    lphi_0[{{ i }}] = {{ constraints_0.lphi_0[i] }};
        {%- endif %}
        {%- if constraints_0.uphi_0[i] != 0 %}
    uphi_0[{{ i }}] = {{ constraints_0.uphi_0[i] }};
        {%- endif %}
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "lphi", lphi_0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "uphi", uphi_0);
    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0,
                                  "nl_constr_phi_o_r_fun", &capsule->phi_0_constraint_fun);
    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0,
                                  "nl_constr_phi_o_r_fun_phi_jac_ux_z_phi_hess_r_jac_ux", &capsule->phi_0_constraint_fun_jac_hess);
    free(luphi_0);
{% endif %}

{% if dims_0.nsh_0 > 0 %}
    // set up soft bounds for nonlinear constraints
    int* idxsh_0 = malloc({{ dims_0.nsh_0 }} * sizeof(int));
    {%- for i in range(end=dims_0.nsh_0) %}
    idxsh_0[{{ i }}] = {{ constraints_0.idxsh_0[i] }};
    {%- endfor %}
    double* lush_0 = calloc(2*{{ dims_0.nsh_0 }}, sizeof(double));
    double* lsh_0 = lush_0;
    double* ush_0 = lush_0 + {{ dims_0.nsh_0 }};
    {%- for i in range(end=dims_0.nsh_0) %}
        {%- if constraints_0.lsh_0[i] != 0 %}
    lsh_0[{{ i }}] = {{ constraints_0.lsh_0[i] }};
        {%- endif %}
        {%- if constraints_0.ush_0[i] != 0 %}
    ush_0[{{ i }}] = {{ constraints_0.ush_0[i] }};
        {%- endif %}
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "idxsh", idxsh_0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "lsh", lsh_0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "ush", ush_0);
    free(idxsh_0);
    free(lush_0);
{%- endif %}

{% if dims_0.nsphi_0 > 0 %}
    // set up soft bounds for convex-over-nonlinear constraints
    int* idxsphi_0 = malloc({{ dims_0.nsphi_0 }} * sizeof(int));
    {%- for i in range(end=dims_0.nsphi_0) %}
    idxsphi_0[{{ i }}] = {{ constraints_0.idxsphi_0[i] }};
    {%- endfor %}
    double* lusphi_0 = calloc(2*{{ dims_0.nsphi_0 }}, sizeof(double));
    double* lsphi_0 = lusphi_0;
    double* usphi_0 = lusphi_0 + {{ dims_0.nsphi_0 }};
    {%- for i in range(end=dims_0.nsphi_0) %}
        {%- if constraints_0.lsphi_0[i] != 0 %}
    lsphi_0[{{ i }}] = {{ constraints_0.lsphi_0[i] }};
        {%- endif %}
        {%- if constraints_0.usphi_0[i] != 0 %}
    usphi_0[{{ i }}] = {{ constraints_0.usphi_0[i] }};
        {%- endif %}
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "idxsphi", idxsphi_0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "lsphi", lsphi_0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "usphi", usphi_0);
    free(idxsphi_0);
    free(lusphi_0);
{%- endif %}

    /* Path related delarations */
    int i_fun;

    // cost
    double* yref;
    double* Vx;
    double* Vu;
    double* Vz;
    double* W;

    // bounds on u
    int* idxbu;
    double* lubu;
    double* lbu;
    double* ubu;

    // bounds on x
    double* lubx;
    double* lbx;
    double* ubx;
    int* idxbx;

    // general linear constraints
    double* D;
    double* C;
    double* lug;
    double* lg;
    double* ug;

    // nonlinear constraints
    double* luh;
    double* lh;
    double* uh;
    double* luphi;
    double* lphi;
    double* uphi;

    // general slack related
    double* zlumem;
    double* Zl;
    double* Zu;
    double* zl;
    double* zu;

    // specific slack types
    int* idxsbx;
    double* lusbx;
    double* lsbx;
    double* usbx;

    int* idxsbu;
    double* lusbu;
    double* lsbu;
    double* usbu;

    int* idxsg;
    double* lusg;
    double* lsg;
    double* usg;

    int* idxsh;
    double* lush;
    double* lsh;
    double* ush;

    int* idxsphi;
    double* lusphi;
    double* lsphi;
    double* usphi;

{%- for jj in range(end=n_phases) %}{# phases loop !#}

    /*********************
     *  Phase {{ jj }}
     * *******************/
    /**** Dynamics ****/
    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        i_fun = i - {{ start_idx[jj] }};
    {%- if mocp_opts.integrator_type[jj] == "ERK" %}
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "expl_vde_forw", &capsule->expl_vde_forw_{{ jj }}[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "expl_vde_adj", &capsule->expl_vde_adj_{{ jj }}[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "expl_ode_fun", &capsule->expl_ode_fun_{{ jj }}[i_fun]);
        {%- if solver_options.hessian_approx == "EXACT" %}
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "expl_ode_hess", &capsule->expl_ode_hess_{{ jj }}[i_fun]);
        {%- endif %}
    {%- elif mocp_opts.integrator_type[jj] == "IRK" %}
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "impl_dae_fun", &capsule->impl_dae_fun_{{ jj }}[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i,
                                   "impl_dae_fun_jac_x_xdot_z", &capsule->impl_dae_fun_jac_x_xdot_z_{{ jj }}[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i,
                                   "impl_dae_jac_x_xdot_u", &capsule->impl_dae_jac_x_xdot_u_z_{{ jj }}[i_fun]);
        {%- if solver_options.hessian_approx == "EXACT" %}
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "impl_dae_hess", &capsule->impl_dae_hess_{{ jj }}[i_fun]);
        {%- endif %}
    {%- elif mocp_opts.integrator_type[jj] == "LIFTED_IRK" %}
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "impl_dae_fun", &capsule->impl_dae_fun_{{ jj }}[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i,
                                   "impl_dae_fun_jac_x_xdot_u", &capsule->impl_dae_fun_jac_x_xdot_u_{{ jj }}[i_fun]);
    {%- elif mocp_opts.integrator_type[jj] == "GNSF" %}
        {% if model[jj].gnsf_purely_linear != 1 %}
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "phi_fun", &capsule->gnsf_phi_fun_{{ jj }}[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "phi_fun_jac_y", &capsule->gnsf_phi_fun_jac_y_{{ jj }}[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "phi_jac_y_uhat", &capsule->gnsf_phi_jac_y_uhat_{{ jj }}[i_fun]);
            {% if model[jj].gnsf_nontrivial_f_LO == 1 %}
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "f_lo_jac_x1_x1dot_u_z",
                                   &capsule->gnsf_f_lo_jac_x1_x1dot_u_z_{{ jj }}[i_fun]);
            {%- endif %}
        {%- endif %}
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "gnsf_get_matrices_fun",
                                   &capsule->gnsf_get_matrices_fun_{{ jj }}[i_fun]);
    {%- elif mocp_opts.integrator_type[jj] == "DISCRETE" %}
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "disc_dyn_fun", &capsule->discr_dyn_phi_fun_{{ jj }}[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "disc_dyn_fun_jac",
                                   &capsule->discr_dyn_phi_fun_jac_ut_xt_{{ jj }}[i_fun]);
        {%- if solver_options.hessian_approx == "EXACT" %}
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "disc_dyn_fun_jac_hess",
                                   &capsule->discr_dyn_phi_fun_jac_ut_xt_hess_{{ jj }}[i_fun]);
        {%- endif %}
    {%- endif %}
    }


{%- if cost[jj].cost_type == "NONLINEAR_LS" %}
    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        i_fun = i - {{ cost_start_idx[jj] }};

        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun", &capsule->cost_y_fun_{{ jj }}[i_fun]);
        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun_jac", &capsule->cost_y_fun_jac_ut_xt_{{ jj }}[i_fun]);
        {%- if solver_options.hessian_approx == "EXACT" %}
        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nls_y_hess", &capsule->cost_y_hess_{{ jj }}[i_fun]);
        {%- endif %}
    }
{%- elif cost[jj].cost_type == "CONVEX_OVER_NONLINEAR" %}
    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        i_fun = i - {{ cost_start_idx[jj] }};
        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "conl_cost_fun", &capsule->conl_cost_fun_{{ jj }}[i_fun]);
        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "conl_cost_fun_jac_hess", &capsule->conl_cost_fun_jac_hess_{{ jj }}[i_fun]);
    }
{%- elif cost[jj].cost_type == "EXTERNAL" %}
    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        i_fun = i - {{ cost_start_idx[jj] }};
        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "ext_cost_fun", &capsule->ext_cost_fun_{{ jj }}[i_fun]);
        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "ext_cost_fun_jac", &capsule->ext_cost_fun_jac_{{ jj }}[i_fun]);
        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "ext_cost_fun_jac_hess", &capsule->ext_cost_fun_jac_hess_{{ jj }}[i_fun]);
    }
{%- endif %}



{%- if mocp_opts.cost_discretization[jj] == "INTEGRATOR" %}
    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        i_fun = i - {{ cost_start_idx[jj] }};
  {%- if cost[jj].cost_type == "NONLINEAR_LS" %}
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun_jac", &capsule->cost_y_fun_jac_ut_xt_{{ jj }}[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun", &capsule->cost_y_fun_{{ jj }}[i_fun]);
  {%- elif cost[jj].cost_type == "CONVEX_OVER_NONLINEAR" %}
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "conl_cost_fun", &capsule->conl_cost_fun_{{ jj }}[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "conl_cost_fun_jac_hess", &capsule->conl_cost_fun_jac_hess_{{ jj }}[i_fun]);
  {%- endif %}
    }
{%- endif %}

    /**** Cost phase {{ jj }} ****/

{%- if phases_dims[jj].ny != 0 %}
    yref = calloc({{ phases_dims[jj].ny }}, sizeof(double));
    // change only the non-zero elements:
    {%- for j in range(end=phases_dims[jj].ny) %}
        {%- if cost[jj].yref[j] != 0 %}
    yref[{{ j }}] = {{ cost[jj].yref[j] }};
        {%- endif %}
    {%- endfor %}

    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "yref", yref);
    }
    free(yref);
  {%- if cost[jj].cost_type == "NONLINEAR_LS" or cost[jj].cost_type == "LINEAR_LS" %}
    W = calloc({{ phases_dims[jj].ny }}*{{ phases_dims[jj].ny }}, sizeof(double));
    // change only the non-zero elements:
    {%- for j in range(end=phases_dims[jj].ny) %}
        {%- for k in range(end=phases_dims[jj].ny) %}
            {%- if cost[jj].W[j][k] != 0 %}
    W[{{ j }}+({{ phases_dims[jj].ny }}) * {{ k }}] = {{ cost[jj].W[j][k] }};
            {%- endif %}
        {%- endfor %}
    {%- endfor %}

    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "W", W);
    }
    free(W);
  {%- endif %}

  {%- if cost[jj].cost_type == "LINEAR_LS" %}
    Vx = calloc({{ phases_dims[jj].ny }}*{{ phases_dims[jj].nx }}, sizeof(double));
    // change only the non-zero elements:
    {%- for j in range(end=phases_dims[jj].ny) %}
        {%- for k in range(end=phases_dims[jj].nx) %}
            {%- if cost[jj].Vx[j][k] != 0 %}
    Vx[{{ j }}+({{ phases_dims[jj].ny }}) * {{ k }}] = {{ cost[jj].Vx[j][k] }};
            {%- endif %}
        {%- endfor %}
    {%- endfor %}
    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Vx", Vx);
    }
    free(Vx);

    {% if phases_dims[jj].nu > 0 %}
    Vu = calloc({{ phases_dims[jj].ny }}*{{ phases_dims[jj].nu }}, sizeof(double));
    // change only the non-zero elements:
    {%- for j in range(end=phases_dims[jj].ny) %}
        {%- for k in range(end=phases_dims[jj].nu) %}
            {%- if cost[jj].Vu[j][k] != 0 %}
    Vu[{{ j }}+({{ phases_dims[jj].ny }}) * {{ k }}] = {{ cost[jj].Vu[j][k] }};
            {%- endif %}
        {%- endfor %}
    {%- endfor %}

    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Vu", Vu);
    }
    free(Vu);
    {%- endif %}

    {%- if phases_dims[jj].nz > 0 %}
    Vz = calloc({{ phases_dims[jj].ny }}*{{ phases_dims[jj].nz }}, sizeof(double));
    // change only the non-zero elements:
    {%- for j in range(end=phases_dims[jj].ny) %}
        {%- for k in range(end=phases_dims[jj].nz) %}
            {%- if cost[jj].Vz[j][k] != 0 %}
    Vz[{{ j }}+({{ phases_dims[jj].ny }}) * {{ k }}] = {{ cost[jj].Vz[j][k] }};
            {%- endif %}
        {%- endfor %}
    {%- endfor %}

    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Vz", Vz);
    }
    free(Vz);
    {%- endif %}
  {%- endif %}{# LINEAR LS #}
{%- endif %}{# ny != 0 #}


{%- if phases_dims[jj].ns > 0 %}
    zlumem = calloc(4*{{ phases_dims[jj].ns }}, sizeof(double));
    Zl = zlumem+{{ phases_dims[jj].ns }}*0;
    Zu = zlumem+{{ phases_dims[jj].ns }}*1;
    zl = zlumem+{{ phases_dims[jj].ns }}*2;
    zu = zlumem+{{ phases_dims[jj].ns }}*3;
    // change only the non-zero elements:
    {%- for j in range(end=phases_dims[jj].ns) %}
        {%- if cost[jj].Zl[j] != 0 %}
    Zl[{{ j }}] = {{ cost[jj].Zl[j] }};
        {%- endif %}
    {%- endfor %}

    {%- for j in range(end=phases_dims[jj].ns) %}
        {%- if cost[jj].Zu[j] != 0 %}
    Zu[{{ j }}] = {{ cost[jj].Zu[j] }};
        {%- endif %}
    {%- endfor %}

    {%- for j in range(end=phases_dims[jj].ns) %}
        {%- if cost[jj].zl[j] != 0 %}
    zl[{{ j }}] = {{ cost[jj].zl[j] }};
        {%- endif %}
    {%- endfor %}

    {%- for j in range(end=phases_dims[jj].ns) %}
        {%- if cost[jj].zu[j] != 0 %}
    zu[{{ j }}] = {{ cost[jj].zu[j] }};
        {%- endif %}
    {%- endfor %}

    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Zl", Zl);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Zu", Zu);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "zl", zl);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "zu", zu);
    }
    free(zlumem);
{%- endif %}


    /**** Constraints phase {{ jj }} ****/

    /* constraints that are the same for initial and intermediate */
{%- if phases_dims[jj].nbu > 0 %}
    // u
    idxbu = malloc({{ phases_dims[jj].nbu }} * sizeof(int));
    {%- for i in range(end=phases_dims[jj].nbu) %}
    idxbu[{{ i }}] = {{ constraints[jj].idxbu[i] }};
    {%- endfor %}
    lubu = calloc(2*{{ phases_dims[jj].nbu }}, sizeof(double));
    lbu = lubu;
    ubu = lubu + {{ phases_dims[jj].nbu }};
    {%- for i in range(end=phases_dims[jj].nbu) %}
        {%- if constraints[jj].lbu[i] != 0 %}
    lbu[{{ i }}] = {{ constraints[jj].lbu[i] }};
        {%- endif %}
        {%- if constraints[jj].ubu[i] != 0 %}
    ubu[{{ i }}] = {{ constraints[jj].ubu[i] }};
        {%- endif %}
    {%- endfor %}

    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "idxbu", idxbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "lbu", lbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "ubu", ubu);
    }
    free(idxbu);
    free(lubu);
{%- endif %}


{% if phases_dims[jj].ng > 0 %}
    // set up general constraints for stage 0 to N-1
    D = calloc({{ phases_dims[jj].ng }}*{{ phases_dims[jj].nu }}, sizeof(double));
    C = calloc({{ phases_dims[jj].ng }}*{{ phases_dims[jj].nx }}, sizeof(double));
    lug = calloc(2*{{ phases_dims[jj].ng }}, sizeof(double));
    lg = lug;
    ug = lug + {{ phases_dims[jj].ng }};

    {%- for j in range(end=phases_dims[jj].ng) -%}
        {% for k in range(end=phases_dims[jj].nu) %}
            {%- if constraints[jj].D[j][k] != 0 %}
    D[{{ j }}+{{ phases_dims[jj].ng }} * {{ k }}] = {{ constraints[jj].D[j][k] }};
            {%- endif %}
        {%- endfor %}
    {%- endfor %}

    {%- for j in range(end=phases_dims[jj].ng) -%}
        {% for k in range(end=phases_dims[jj].nx) %}
            {%- if constraints[jj].C[j][k] != 0 %}
    C[{{ j }}+{{ phases_dims[jj].ng }} * {{ k }}] = {{ constraints[jj].C[j][k] }};
            {%- endif %}
        {%- endfor %}
    {%- endfor %}

    {%- for i in range(end=phases_dims[jj].ng) %}
        {%- if constraints[jj].lg[i] != 0 %}
    lg[{{ i }}] = {{ constraints[jj].lg[i] }};
        {%- endif %}
    {%- endfor %}

    {%- for i in range(end=phases_dims[jj].ng) %}
        {%- if constraints[jj].ug[i] != 0 %}
    ug[{{ i }}] = {{ constraints[jj].ug[i] }};
        {%- endif %}
    {%- endfor %}

    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "D", D);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "C", C);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "lg", lg);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "ug", ug);
    }
    free(D);
    free(C);
    free(lug);
{%- endif %}


{%- if phases_dims[jj].nsbu > 0 %}
    // set up soft bounds for u
    idxsbu = malloc({{ phases_dims[jj].nsbu }} * sizeof(int));
    {%- for i in range(end=phases_dims[jj].nsbu) %}
    idxsbu[{{ i }}] = {{ constraints[jj].idxsbu[i] }};
    {%- endfor %}
    lusbu = calloc(2*{{ phases_dims[jj].nsbu }}, sizeof(double));
    lsbu = lusbu;
    usbu = lusbu + {{ phases_dims[jj].nsbu }};
    {%- for i in range(end=phases_dims[jj].nsbu) %}
        {%- if constraints[jj].lsbu[i] != 0 %}
    lsbu[{{ i }}] = {{ constraints[jj].lsbu[i] }};
        {%- endif %}
        {%- if constraints[jj].usbu[i] != 0 %}
    usbu[{{ i }}] = {{ constraints[jj].usbu[i] }};
        {%- endif %}
    {%- endfor %}
    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "idxsbu", idxsbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "lsbu", lsbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "usbu", usbu);
    }
    free(idxsbu);
    free(lusbu);
{%- endif %}

{% if phases_dims[jj].nsg > 0 %}
    // set up soft bounds for general linear constraints
    idxsg = malloc({{ phases_dims[jj].nsg }} * sizeof(int));
    {%- for i in range(end=phases_dims[jj].nsg) %}
    idxsg[{{ i }}] = {{ constraints[jj].idxsg[i] }};
    {%- endfor %}
    lusg = calloc(2*{{ phases_dims[jj].nsg }}, sizeof(double));
    lsg = lusg;
    usg = lusg + {{ phases_dims[jj].nsg }};
    {%- for i in range(end=phases_dims[jj].nsg) %}
        {%- if constraints[jj].lsg[i] != 0 %}
    lsg[{{ i }}] = {{ constraints[jj].lsg[i] }};
        {%- endif %}
        {%- if constraints[jj].usg[i] != 0 %}
    usg[{{ i }}] = {{ constraints[jj].usg[i] }};
        {%- endif %}
    {%- endfor %}

    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "idxsg", idxsg);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "lsg", lsg);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "usg", usg);
    }
    free(idxsg);
    free(lusg);
{%- endif %}

    /* Path constraints */
{% if phases_dims[jj].nbx > 0 %}
    // x
    idxbx = malloc({{ phases_dims[jj].nbx }} * sizeof(int));
    {%- for i in range(end=phases_dims[jj].nbx) %}
    idxbx[{{ i }}] = {{ constraints[jj].idxbx[i] }};
    {%- endfor %}
    lubx = calloc(2*{{ phases_dims[jj].nbx }}, sizeof(double));
    lbx = lubx;
    ubx = lubx + {{ phases_dims[jj].nbx }};
    {%- for i in range(end=phases_dims[jj].nbx) %}
        {%- if constraints[jj].lbx[i] != 0 %}
    lbx[{{ i }}] = {{ constraints[jj].lbx[i] }};
        {%- endif %}
        {%- if constraints[jj].ubx[i] != 0 %}
    ubx[{{ i }}] = {{ constraints[jj].ubx[i] }};
        {%- endif %}
    {%- endfor %}

    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "idxbx", idxbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "lbx", lbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "ubx", ubx);
    }
    free(idxbx);
    free(lubx);
{%- endif %}

{% if phases_dims[jj].nh > 0 %}
    // set up nonlinear constraints for stage 1 to N-1
    luh = calloc(2*{{ phases_dims[jj].nh }}, sizeof(double));
    lh = luh;
    uh = luh + {{ phases_dims[jj].nh }};

    {%- for i in range(end=phases_dims[jj].nh) %}
        {%- if constraints[jj].lh[i] != 0 %}
    lh[{{ i }}] = {{ constraints[jj].lh[i] }};
        {%- endif %}
    {%- endfor %}

    {%- for i in range(end=phases_dims[jj].nh) %}
        {%- if constraints[jj].uh[i] != 0 %}
    uh[{{ i }}] = {{ constraints[jj].uh[i] }};
        {%- endif %}
    {%- endfor %}

    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        i_fun = i - {{ cost_start_idx[jj] }};
        ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nl_constr_h_fun_jac",
                                      &capsule->nl_constr_h_fun_jac_{{ jj }}[i_fun]);
        ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nl_constr_h_fun",
                                      &capsule->nl_constr_h_fun_{{ jj }}[i_fun]);
        {% if solver_options.hessian_approx == "EXACT" %}
        ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i,
                                      "nl_constr_h_fun_jac_hess", &capsule->nl_constr_h_fun_jac_hess_{{ jj }}[i_fun]);
        {% endif %}
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "lh", lh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "uh", uh);
        {% if solver_options.with_solution_sens_wrt_params %}
        ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i,
                                      "nl_constr_h_jac_p_hess_xu_p", &capsule->nl_constr_h_jac_p_hess_xu_p_{{ jj }}[i_fun]);
        {% endif %}
        {% if solver_options.with_value_sens_wrt_params %}
        ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i,
                                      "nl_constr_h_adj_p", &capsule->nl_constr_h_adj_p_{{ jj }}[i_fun]);
        {% endif %}
    }
    free(luh);
{%- endif %}

{% if phases_dims[jj].nphi > 0 and constraints[jj].constr_type == "BGP" %}
    // set up convex-over-nonlinear constraints for stage 1 to N-1
    luphi = calloc(2*{{ phases_dims[jj].nphi }}, sizeof(double));
    lphi = luphi;
    uphi = luphi + {{ phases_dims[jj].nphi }};
    {%- for i in range(end=phases_dims[jj].nphi) %}
        {%- if constraints[jj].lphi[i] != 0 %}
    lphi[{{ i }}] = {{ constraints[jj].lphi[i] }};
        {%- endif %}
    {%- endfor %}

    {%- for i in range(end=phases_dims[jj].nphi) %}
        {%- if constraints[jj].uphi[i] != 0 %}
    uphi[{{ i }}] = {{ constraints[jj].uphi[i] }};
        {%- endif %}
    {%- endfor %}

    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        i_fun = i - {{ cost_start_idx[jj] }};
        ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nl_constr_phi_o_r_fun", &capsule->phi_constraint_fun[i_fun]);
        ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i,
                                      "nl_constr_phi_o_r_fun_phi_jac_ux_z_phi_hess_r_jac_ux", &capsule->phi_constraint_fun_jac_hess[i_fun]);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "lphi", lphi);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "uphi", uphi);
    }
    free(luphi);
{%- endif %}


{%- if phases_dims[jj].nsbx > 0 %}
{# TODO: introduce nsbx0 #}
    // ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "idxsbx", idxsbx);
    // ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "lsbx", lsbx);
    // ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "usbx", usbx);

    // soft bounds on x
    idxsbx = malloc({{ phases_dims[jj].nsbx }} * sizeof(int));
    {%- for i in range(end=phases_dims[jj].nsbx) %}
    idxsbx[{{ i }}] = {{ constraints[jj].idxsbx[i] }};
    {%- endfor %}

    lusbx = calloc(2*{{ phases_dims[jj].nsbx }}, sizeof(double));
    lsbx = lusbx;
    usbx = lusbx + {{ phases_dims[jj].nsbx }};
    {%- for i in range(end=phases_dims[jj].nsbx) %}
        {%- if constraints[jj].lsbx[i] != 0 %}
    lsbx[{{ i }}] = {{ constraints[jj].lsbx[i] }};
        {%- endif %}
        {%- if constraints[jj].usbx[i] != 0 %}
    usbx[{{ i }}] = {{ constraints[jj].usbx[i] }};
        {%- endif %}
    {%- endfor %}

    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "idxsbx", idxsbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "lsbx", lsbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "usbx", usbx);
    }
    free(idxsbx);
    free(lusbx);
{%- endif %}


{% if phases_dims[jj].nsh > 0 %}
    // set up soft bounds for nonlinear constraints
    idxsh = malloc({{ phases_dims[jj].nsh }} * sizeof(int));
    {%- for i in range(end=phases_dims[jj].nsh) %}
    idxsh[{{ i }}] = {{ constraints[jj].idxsh[i] }};
    {%- endfor %}
    lush = calloc(2*{{ phases_dims[jj].nsh }}, sizeof(double));
    lsh = lush;
    ush = lush + {{ phases_dims[jj].nsh }};
    {%- for i in range(end=phases_dims[jj].nsh) %}
        {%- if constraints[jj].lsh[i] != 0 %}
    lsh[{{ i }}] = {{ constraints[jj].lsh[i] }};
        {%- endif %}
        {%- if constraints[jj].ush[i] != 0 %}
    ush[{{ i }}] = {{ constraints[jj].ush[i] }};
        {%- endif %}
    {%- endfor %}

    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "idxsh", idxsh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "lsh", lsh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "ush", ush);
    }
    free(idxsh);
    free(lush);
{%- endif %}

{% if phases_dims[jj].nsphi > 0 %}
    // set up soft bounds for convex-over-nonlinear constraints
    idxsphi = malloc({{ phases_dims[jj].nsphi }} * sizeof(int));
    {%- for i in range(end=phases_dims[jj].nsphi) %}
    idxsphi[{{ i }}] = {{ constraints[jj].idxsphi[i] }};
    {%- endfor %}
    lusphi = calloc(2*{{ phases_dims[jj].nsphi }}, sizeof(double));
    lsphi = lusphi;
    usphi = lusphi + {{ phases_dims[jj].nsphi }};
    {%- for i in range(end=phases_dims[jj].nsphi) %}
        {%- if constraints[jj].lsphi[i] != 0 %}
    lsphi[{{ i }}] = {{ constraints[jj].lsphi[i] }};
        {%- endif %}
        {%- if constraints[jj].usphi[i] != 0 %}
    usphi[{{ i }}] = {{ constraints[jj].usphi[i] }};
        {%- endif %}
    {%- endfor %}

    for (int i = {{ cost_start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "idxsphi", idxsphi);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "lsphi", lsphi);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i, "usphi", usphi);
    }
    free(idxsphi);
    free(lusphi);
{%- endif %}

{%- endfor %}{# phases loop !#}



    // TERMINAL node
{%- if dims_e.ny_e != 0 %}
    double* yref_e = calloc({{ dims_e.ny_e }}, sizeof(double));
    // change only the non-zero elements:
    {%- for j in range(end=dims_e.ny_e) %}
        {%- if cost_e.yref_e[j] != 0 %}
    yref_e[{{ j }}] = {{ cost_e.yref_e[j] }};
        {%- endif %}
    {%- endfor %}
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "yref", yref_e);
    free(yref_e);

  {%- if cost_e.cost_type_e == "NONLINEAR_LS" or cost_e.cost_type_e == "LINEAR_LS" %}

    double* W_e = calloc({{ dims_e.ny_e }}*{{ dims_e.ny_e }}, sizeof(double));
    // change only the non-zero elements:
    {%- for j in range(end=dims_e.ny_e) %}
        {%- for k in range(end=dims_e.ny_e) %}
            {%- if cost_e.W_e[j][k] != 0 %}
    W_e[{{ j }}+({{ dims_e.ny_e }}) * {{ k }}] = {{ cost_e.W_e[j][k] }};
            {%- endif %}
        {%- endfor %}
    {%- endfor %}
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "W", W_e);
    free(W_e);
  {%- endif %}


  {%- if cost_e.cost_type_e == "LINEAR_LS" %}
    double* Vx_e = calloc({{ dims_e.ny_e }}*{{ dims_e.nx }}, sizeof(double));
    // change only the non-zero elements:
    {%- for j in range(end=dims_e.ny_e) %}
        {%- for k in range(end=dims_e.nx) %}
            {%- if cost_e.Vx_e[j][k] != 0 %}
    Vx_e[{{ j }}+({{ dims_e.ny_e }}) * {{ k }}] = {{ cost_e.Vx_e[j][k] }};
            {%- endif %}
        {%- endfor %}
    {%- endfor %}
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "Vx", Vx_e);
    free(Vx_e);
  {%- endif %}
{%- endif %}{# ny_e != 0 #}


{%- if cost_e.cost_type_e == "NONLINEAR_LS" %}
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "nls_y_fun", &capsule->cost_y_e_fun);
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "nls_y_fun_jac", &capsule->cost_y_e_fun_jac_ut_xt);
    {%- if solver_options.hessian_approx == "EXACT" %}
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "nls_y_hess", &capsule->cost_y_e_hess);
    {%- endif %}

{%- elif cost_e.cost_type_e == "CONVEX_OVER_NONLINEAR" %}
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "conl_cost_fun", &capsule->conl_cost_e_fun);
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "conl_cost_fun_jac_hess", &capsule->conl_cost_e_fun_jac_hess);

{%- elif cost_e.cost_type_e == "EXTERNAL" %}
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "ext_cost_fun", &capsule->ext_cost_e_fun);
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "ext_cost_fun_jac", &capsule->ext_cost_e_fun_jac);
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "ext_cost_fun_jac_hess", &capsule->ext_cost_e_fun_jac_hess);
{%- endif %}



{% if dims_e.ns_e > 0 %}
    // slacks terminal
    double* zluemem = calloc(4*{{ dims_e.ns_e }}, sizeof(double));
    double* Zl_e = zluemem+{{ dims_e.ns_e }}*0;
    double* Zu_e = zluemem+{{ dims_e.ns_e }}*1;
    double* zl_e = zluemem+{{ dims_e.ns_e }}*2;
    double* zu_e = zluemem+{{ dims_e.ns_e }}*3;

    // change only the non-zero elements:
    {%- for j in range(end=dims_e.ns_e) %}
        {%- if cost_e.Zl_e[j] != 0 %}
    Zl_e[{{ j }}] = {{ cost_e.Zl_e[j] }};
        {%- endif %}
    {%- endfor %}

    {%- for j in range(end=dims_e.ns_e) %}
        {%- if cost_e.Zu_e[j] != 0 %}
    Zu_e[{{ j }}] = {{ cost_e.Zu_e[j] }};
        {%- endif %}
    {%- endfor %}

    {%- for j in range(end=dims_e.ns_e) %}
        {%- if cost_e.zl_e[j] != 0 %}
    zl_e[{{ j }}] = {{ cost_e.zl_e[j] }};
        {%- endif %}
    {%- endfor %}

    {%- for j in range(end=dims_e.ns_e) %}
        {%- if cost_e.zu_e[j] != 0 %}
    zu_e[{{ j }}] = {{ cost_e.zu_e[j] }};
        {%- endif %}
    {%- endfor %}

    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "Zl", Zl_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "Zu", Zu_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "zl", zl_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "zu", zu_e);
    free(zluemem);
{%- endif %}





    /* terminal constraints */
{% if dims_e.nbx_e > 0 %}
    // set up bounds for last stage
    // x
    int* idxbx_e = malloc({{ dims_e.nbx_e }} * sizeof(int));
    {%- for i in range(end=dims_e.nbx_e) %}
    idxbx_e[{{ i }}] = {{ constraints_e.idxbx_e[i] }};
    {%- endfor %}
    double* lubx_e = calloc(2*{{ dims_e.nbx_e }}, sizeof(double));
    double* lbx_e = lubx_e;
    double* ubx_e = lubx_e + {{ dims_e.nbx_e }};
    {%- for i in range(end=dims_e.nbx_e) %}
        {%- if constraints_e.lbx_e[i] != 0 %}
    lbx_e[{{ i }}] = {{ constraints_e.lbx_e[i] }};
        {%- endif %}
        {%- if constraints_e.ubx_e[i] != 0 %}
    ubx_e[{{ i }}] = {{ constraints_e.ubx_e[i] }};
        {%- endif %}
    {%- endfor %}
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "idxbx", idxbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "lbx", lbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "ubx", ubx_e);
    free(idxbx_e);
    free(lubx_e);
{%- endif %}

{% if dims_e.nsg_e > 0 %}
    // set up soft bounds for general linear constraints
    int* idxsg_e = calloc({{ dims_e.nsg_e }}, sizeof(int));
    {%- for i in range(end=dims_e.nsg_e) %}
    idxsg_e[{{ i }}] = {{ constraints_e.idxsg_e[i] }};
    {%- endfor %}
    double* lusg_e = calloc(2*{{ dims_e.nsg_e }}, sizeof(double));
    double* lsg_e = lusg_e;
    double* usg_e = lusg_e + {{ dims_e.nsg_e }};
    {%- for i in range(end=dims_e.nsg_e) %}
        {%- if constraints_e.lsg_e[i] != 0 %}
    lsg_e[{{ i }}] = {{ constraints_e.lsg_e[i] }};
        {%- endif %}
        {%- if constraints_e.usg_e[i] != 0 %}
    usg_e[{{ i }}] = {{ constraints_e.usg_e[i] }};
        {%- endif %}
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "idxsg", idxsg_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "lsg", lsg_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "usg", usg_e);
    free(idxsg_e);
    free(lusg_e);
{%- endif %}

{% if dims_e.nsh_e > 0 %}
    // set up soft bounds for nonlinear constraints
    int* idxsh_e = malloc({{ dims_e.nsh_e }} * sizeof(int));
    {%- for i in range(end=dims_e.nsh_e) %}
    idxsh_e[{{ i }}] = {{ constraints_e.idxsh_e[i] }};
    {%- endfor %}
    double* lush_e = calloc(2*{{ dims_e.nsh_e }}, sizeof(double));
    double* lsh_e = lush_e;
    double* ush_e = lush_e + {{ dims_e.nsh_e }};
    {%- for i in range(end=dims_e.nsh_e) %}
        {%- if constraints_e.lsh_e[i] != 0 %}
    lsh_e[{{ i }}] = {{ constraints_e.lsh_e[i] }};
        {%- endif %}
        {%- if constraints_e.ush_e[i] != 0 %}
    ush_e[{{ i }}] = {{ constraints_e.ush_e[i] }};
        {%- endif %}
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "idxsh", idxsh_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "lsh", lsh_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "ush", ush_e);
    free(idxsh_e);
    free(lush_e);
{%- endif %}

{% if dims_e.nsphi_e > 0 %}
    // set up soft bounds for convex-over-nonlinear constraints
    int* idxsphi_e = malloc({{ dims_e.nsphi_e }} * sizeof(int));
    {%- for i in range(end=dims_e.nsphi_e) %}
    idxsphi_e[{{ i }}] = {{ constraints_e.idxsphi_e[i] }};
    {%- endfor %}
    double* lusphi_e = calloc(2*{{ dims_e.nsphi_e }}, sizeof(double));
    double* lsphi_e = lusphi_e;
    double* usphi_e = lusphi_e + {{ dims_e.nsphi_e }};
    {%- for i in range(end=dims_e.nsphi_e) %}
        {%- if constraints_e.lsphi_e[i] != 0 %}
    lsphi_e[{{ i }}] = {{ constraints_e.lsphi_e[i] }};
        {%- endif %}
        {%- if constraints_e.usphi_e[i] != 0 %}
    usphi_e[{{ i }}] = {{ constraints_e.usphi_e[i] }};
        {%- endif %}
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "idxsphi", idxsphi_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "lsphi", lsphi_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "usphi", usphi_e);
    free(idxsphi_e);
    free(lusphi_e);
{%- endif %}

{% if dims_e.nsbx_e > 0 %}
    // soft bounds on x
    int* idxsbx_e = malloc({{ dims_e.nsbx_e }} * sizeof(int));
    {%- for i in range(end=dims_e.nsbx_e) %}
    idxsbx_e[{{ i }}] = {{ constraints_e.idxsbx_e[i] }};
    {%- endfor %}
    double* lusbx_e = calloc(2*{{ dims_e.nsbx_e }}, sizeof(double));
    double* lsbx_e = lusbx_e;
    double* usbx_e = lusbx_e + {{ dims_e.nsbx_e }};
    {%- for i in range(end=dims_e.nsbx_e) %}
        {%- if constraints_e.lsbx_e[i] != 0 %}
    lsbx_e[{{ i }}] = {{ constraints_e.lsbx_e[i] }};
        {%- endif %}
        {%- if constraints_e.usbx_e[i] != 0 %}
    usbx_e[{{ i }}] = {{ constraints_e.usbx_e[i] }};
        {%- endif %}
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "idxsbx", idxsbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "lsbx", lsbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "usbx", usbx_e);
    free(idxsbx_e);
    free(lusbx_e);
{% endif %}

{% if dims_e.ng_e > 0 %}
    // set up general constraints for last stage
    double* C_e = calloc({{ dims_e.ng_e }}*{{ dims_e.nx }}, sizeof(double));
    double* lug_e = calloc(2*{{ dims_e.ng_e }}, sizeof(double));
    double* lg_e = lug_e;
    double* ug_e = lug_e + {{ dims_e.ng_e }};

    {%- for j in range(end=dims_e.ng_e) %}
        {%- for k in range(end=dims_e.nx) %}
            {%- if constraints_e.C_e[j][k] != 0 %}
    C_e[{{ j }}+{{ dims_e.ng_e }} * {{ k }}] = {{ constraints_e.C_e[j][k] }};
            {%- endif %}
        {%- endfor %}
    {%- endfor %}

    {%- for i in range(end=dims_e.ng_e) %}
        {%- if constraints_e.lg_e[i] != 0 %}
    lg_e[{{ i }}] = {{ constraints_e.lg_e[i] }};
        {%- endif %}
        {%- if constraints_e.ug_e[i] != 0 %}
    ug_e[{{ i }}] = {{ constraints_e.ug_e[i] }};
        {%- endif %}
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "C", C_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "lg", lg_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "ug", ug_e);
    free(C_e);
    free(lug_e);
{%- endif %}

{% if dims_e.nh_e > 0 %}
    // set up nonlinear constraints for last stage
    double* luh_e = calloc(2*{{ dims_e.nh_e }}, sizeof(double));
    double* lh_e = luh_e;
    double* uh_e = luh_e + {{ dims_e.nh_e }};
    {%- for i in range(end=dims_e.nh_e) %}
        {%- if constraints_e.lh_e[i] != 0 %}
    lh_e[{{ i }}] = {{ constraints_e.lh_e[i] }};
        {%- endif %}
    {%- endfor %}

    {%- for i in range(end=dims_e.nh_e) %}
        {%- if constraints_e.uh_e[i] != 0 %}
    uh_e[{{ i }}] = {{ constraints_e.uh_e[i] }};
        {%- endif %}
    {%- endfor %}

    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "nl_constr_h_fun_jac", &capsule->nl_constr_h_e_fun_jac);
    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "nl_constr_h_fun", &capsule->nl_constr_h_e_fun);
    {% if solver_options.hessian_approx == "EXACT" %}
    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "nl_constr_h_fun_jac_hess",
                                  &capsule->nl_constr_h_e_fun_jac_hess);
    {% endif %}
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "lh", lh_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "uh", uh_e);
    {% if solver_options.with_solution_sens_wrt_params %}
    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "nl_constr_h_jac_p_hess_xu_p",
                                  &capsule->nl_constr_h_e_jac_p_hess_xu_p);
    {% endif %}
    {% if solver_options.with_value_sens_wrt_params %}
    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "nl_constr_h_adj_p",
                                  &capsule->nl_constr_h_e_adj_p);
    {% endif %}
    free(luh_e);
{%- elif dims_e.nphi_e > 0 and constraints_e.constr_type_e == "BGP" %}
    // set up convex-over-nonlinear constraints for last stage
    double* luphi_e = calloc(2*{{ dims_e.nphi_e }}, sizeof(double));
    double* lphi_e = luphi_e;
    double* uphi_e = luphi_e + {{ dims_e.nphi_e }};
    {%- for i in range(end=dims_e.nphi_e) %}
        {%- if constraints_e.lphi_e[i] != 0 %}
    lphi_e[{{ i }}] = {{ constraints_e.lphi_e[i] }};
        {%- endif %}
        {%- if constraints_e.uphi_e[i] != 0 %}
    uphi_e[{{ i }}] = {{ constraints_e.uphi_e[i] }};
        {%- endif %}
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "lphi", lphi_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "uphi", uphi_e);
    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "nl_constr_phi_o_r_fun", &capsule->phi_e_constraint_fun);
    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N,
                                  "nl_constr_phi_o_r_fun_phi_jac_ux_z_phi_hess_r_jac_ux", &capsule->phi_e_constraint_fun_jac_hess);
    free(luphi_e);
{% endif %}
}


/**
 * Internal function for {{ name }}_acados_create: step 6
 */
void {{ name }}_acados_create_set_opts({{ name }}_solver_capsule* capsule)
{
    const int N = capsule->nlp_solver_plan->N;
    ocp_nlp_config* nlp_config = capsule->nlp_config;
    void *nlp_opts = capsule->nlp_opts;

    // declare
    bool tmp_bool;
    int newton_iter_val;
    double newton_tol_val;
    /************************************************
    *  opts
    ************************************************/

{% if solver_options.hessian_approx == "EXACT" %}
    int nlp_solver_exact_hessian = 1;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "exact_hess", &nlp_solver_exact_hessian);

    int exact_hess_dyn = {{ solver_options.exact_hess_dyn }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "exact_hess_dyn", &exact_hess_dyn);

    int exact_hess_cost = {{ solver_options.exact_hess_cost }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "exact_hess_cost", &exact_hess_cost);

    int exact_hess_constr = {{ solver_options.exact_hess_constr }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "exact_hess_constr", &exact_hess_constr);
{%- endif %}

    int fixed_hess = {{ solver_options.fixed_hess }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "fixed_hess", &fixed_hess);

{%- if solver_options.globalization == "FIXED_STEP" %}
    double globalization_fixed_step_length = {{ solver_options.globalization_fixed_step_length }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "globalization_fixed_step_length", &globalization_fixed_step_length);
{% else %}
    double globalization_alpha_min = {{ solver_options.globalization_alpha_min }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "globalization_alpha_min", &globalization_alpha_min);

    double globalization_alpha_reduction = {{ solver_options.globalization_alpha_reduction }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "globalization_alpha_reduction", &globalization_alpha_reduction);
{%- endif %}


{# globalization specific options #}
{%- if solver_options.globalization == "MERIT_BACKTRACKING" %}
    // ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "globalization", "merit_backtracking");

    int globalization_line_search_use_sufficient_descent = {{ solver_options.globalization_line_search_use_sufficient_descent }};
    ocp_nlp_solver_opts_set(nlp_config, capsule->nlp_opts, "globalization_line_search_use_sufficient_descent", &globalization_line_search_use_sufficient_descent);

    int globalization_use_SOC = {{ solver_options.globalization_use_SOC }};
    ocp_nlp_solver_opts_set(nlp_config, capsule->nlp_opts, "globalization_use_SOC", &globalization_use_SOC);

    double globalization_eps_sufficient_descent = {{ solver_options.globalization_eps_sufficient_descent }};
    ocp_nlp_solver_opts_set(nlp_config, capsule->nlp_opts, "globalization_eps_sufficient_descent", &globalization_eps_sufficient_descent);
{%- elif solver_options.globalization == "FUNNEL_L1PEN_LINESEARCH" %}
    // ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "globalization", "funnel_l1pen_linesearch");

    double globalization_funnel_init_increase_factor = {{ solver_options.globalization_funnel_init_increase_factor }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "globalization_funnel_init_increase_factor", &globalization_funnel_init_increase_factor);

    double globalization_funnel_init_upper_bound = {{ solver_options.globalization_funnel_init_upper_bound }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "globalization_funnel_init_upper_bound", &globalization_funnel_init_upper_bound);

    double globalization_funnel_sufficient_decrease_factor = {{ solver_options.globalization_funnel_sufficient_decrease_factor }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "globalization_funnel_sufficient_decrease_factor", &globalization_funnel_sufficient_decrease_factor);

    double globalization_funnel_kappa = {{ solver_options.globalization_funnel_kappa }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "globalization_funnel_kappa", &globalization_funnel_kappa);

    double globalization_funnel_initial_penalty_parameter = {{ solver_options.globalization_funnel_initial_penalty_parameter }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "globalization_funnel_initial_penalty_parameter", &globalization_funnel_initial_penalty_parameter);

    bool globalization_funnel_use_merit_fun_only = {{ solver_options.globalization_funnel_use_merit_fun_only }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "globalization_funnel_use_merit_fun_only", &globalization_funnel_use_merit_fun_only);
{%- endif %}

    int with_solution_sens_wrt_params = {{ solver_options.with_solution_sens_wrt_params }};
    ocp_nlp_solver_opts_set(nlp_config, capsule->nlp_opts, "with_solution_sens_wrt_params", &with_solution_sens_wrt_params);

    int with_value_sens_wrt_params = {{ solver_options.with_value_sens_wrt_params }};
    ocp_nlp_solver_opts_set(nlp_config, capsule->nlp_opts, "with_value_sens_wrt_params", &with_value_sens_wrt_params);

    double solution_sens_qp_t_lam_min = {{ solver_options.solution_sens_qp_t_lam_min }};
    ocp_nlp_solver_opts_set(nlp_config, capsule->nlp_opts, "solution_sens_qp_t_lam_min", &solution_sens_qp_t_lam_min);

    int globalization_full_step_dual = {{ solver_options.globalization_full_step_dual }};
    ocp_nlp_solver_opts_set(nlp_config, capsule->nlp_opts, "globalization_full_step_dual", &globalization_full_step_dual);

    {%- if solver_options.nlp_solver_warm_start_first_qp_from_nlp %}
    int nlp_solver_warm_start_first_qp = {{ solver_options.nlp_solver_warm_start_first_qp }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "warm_start_first_qp", &nlp_solver_warm_start_first_qp);
    {%- endif %}

    {%- if solver_options.nlp_solver_warm_start_first_qp_from_nlp %}
    int nlp_solver_warm_start_first_qp_from_nlp = {{ solver_options.nlp_solver_warm_start_first_qp_from_nlp }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "warm_start_first_qp", &nlp_solver_warm_start_first_qp_from_nlp);
    {%- endif %}

    double levenberg_marquardt = {{ solver_options.levenberg_marquardt }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "levenberg_marquardt", &levenberg_marquardt);

    /* options QP solver */
{%- if solver_options.qp_solver is starting_with("PARTIAL_CONDENSING") %}
    int qp_solver_cond_N;

    {%- if solver_options.qp_solver_cond_N -%}
    const int qp_solver_cond_N_ori = {{ solver_options.qp_solver_cond_N }};
    qp_solver_cond_N = N < qp_solver_cond_N_ori ? N : qp_solver_cond_N_ori; // use the minimum value here
    {%- else %}
    // NOTE: there is no condensing happening here!
    qp_solver_cond_N = N;
    {%- endif %}
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_cond_N", &qp_solver_cond_N);

    {%- if solver_options.qp_solver_cond_block_size -%}
    int* qp_solver_cond_block_size = malloc((qp_solver_cond_N+1) * sizeof(int));
    {%- for i in range(end=solver_options.qp_solver_cond_N+1) %}
    qp_solver_cond_block_size[{{ i }}] = {{ solver_options.qp_solver_cond_block_size[i] }};
    {%- endfor %}
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_cond_block_size", qp_solver_cond_block_size);
    free(qp_solver_cond_block_size);
    {%- endif %}
{%- endif %}

{%- if solver_options.regularize_method == "PROJECT" or solver_options.regularize_method == "MIRROR" or solver_options.regularize_method == "CONVEXIFY" or solver_options.regularize_method == "GERSHGORIN_LEVENBERG_MARQUARDT"%}
    double reg_epsilon = {{ solver_options.reg_epsilon }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "reg_epsilon", &reg_epsilon);
{%- endif %}

{%- if solver_options.regularize_method == "PROJECT" or solver_options.regularize_method == "MIRROR" %}
    double reg_max_cond_block = {{ solver_options.reg_max_cond_block }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "reg_max_cond_block", &reg_max_cond_block);

    double reg_min_epsilon = {{ solver_options.reg_min_epsilon }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "reg_min_epsilon", &reg_min_epsilon);

    bool reg_adaptive_eps = {{ solver_options.reg_adaptive_eps }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "reg_adaptive_eps", &reg_adaptive_eps);
{%- endif %}

    int nlp_solver_ext_qp_res = {{ solver_options.nlp_solver_ext_qp_res }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "ext_qp_res", &nlp_solver_ext_qp_res);

    bool store_iterates = {{ solver_options.store_iterates }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "store_iterates", &store_iterates);

{%- if solver_options.nlp_solver_type == "SQP" or solver_options.nlp_solver_type == "SQP_WITH_FEASIBLE_QP" %}
    int log_primal_step_norm = {{ solver_options.log_primal_step_norm }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "log_primal_step_norm", &log_primal_step_norm);

    int log_dual_step_norm = {{ solver_options.log_dual_step_norm }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "log_dual_step_norm", &log_dual_step_norm);

    double nlp_solver_tol_min_step_norm = {{ solver_options.nlp_solver_tol_min_step_norm }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_min_step_norm", &nlp_solver_tol_min_step_norm);
{%- endif %}

{%- if solver_options.qp_solver is containing("HPIPM") %}
    // set HPIPM mode: should be done before setting other QP solver options
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_hpipm_mode", "{{ solver_options.hpipm_mode }}");

{% if solver_options.qp_solver_mu0 > 0 %}
    double qp_solver_mu0 = {{ solver_options.qp_solver_mu0 }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_mu0", &qp_solver_mu0);
{%- endif %}

    int qp_solver_t0_init = {{ solver_options.qp_solver_t0_init }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_t0_init", &qp_solver_t0_init);
{%- endif %}

{% if solver_options.tau_min > 0 %}
    double tau_min = {{ solver_options.tau_min }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tau_min", &tau_min);
{%- endif %}

{%- if solver_options.nlp_solver_type == "SQP_WITH_FEASIBLE_QP" %}
bool use_constraint_hessian_in_feas_qp = {{ solver_options.use_constraint_hessian_in_feas_qp }};
ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "use_constraint_hessian_in_feas_qp", &use_constraint_hessian_in_feas_qp);

int search_direction_mode = {{ solver_options.search_direction_mode }};
ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "search_direction_mode", &search_direction_mode);

bool allow_direction_mode_switch_to_nominal = {{ solver_options.allow_direction_mode_switch_to_nominal }};
ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "allow_direction_mode_switch_to_nominal", &allow_direction_mode_switch_to_nominal);
{%- endif %}

{% if solver_options.nlp_solver_type == "SQP" or solver_options.nlp_solver_type == "DDP" or solver_options.nlp_solver_type == "SQP_WITH_FEASIBLE_QP" %}
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

    // set options for adaptive Levenberg-Marquardt Update
    bool with_adaptive_levenberg_marquardt = {{ solver_options.with_adaptive_levenberg_marquardt }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "with_adaptive_levenberg_marquardt", &with_adaptive_levenberg_marquardt);

    double adaptive_levenberg_marquardt_lam = {{ solver_options.adaptive_levenberg_marquardt_lam }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "adaptive_levenberg_marquardt_lam", &adaptive_levenberg_marquardt_lam);

    double adaptive_levenberg_marquardt_mu_min = {{ solver_options.adaptive_levenberg_marquardt_mu_min }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "adaptive_levenberg_marquardt_mu_min", &adaptive_levenberg_marquardt_mu_min);

    double adaptive_levenberg_marquardt_mu0 = {{ solver_options.adaptive_levenberg_marquardt_mu0 }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "adaptive_levenberg_marquardt_mu0", &adaptive_levenberg_marquardt_mu0);

    double adaptive_levenberg_marquardt_obj_scalar = {{ solver_options.adaptive_levenberg_marquardt_obj_scalar }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "adaptive_levenberg_marquardt_obj_scalar", &adaptive_levenberg_marquardt_obj_scalar);

    bool eval_residual_at_max_iter = {{ solver_options.eval_residual_at_max_iter }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "eval_residual_at_max_iter", &eval_residual_at_max_iter);

    // QP scaling
    double qpscaling_ub_max_abs_eig = {{ solver_options.qpscaling_ub_max_abs_eig }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qpscaling_ub_max_abs_eig", &qpscaling_ub_max_abs_eig);

    double qpscaling_lb_norm_inf_grad_obj = {{ solver_options.qpscaling_lb_norm_inf_grad_obj }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qpscaling_lb_norm_inf_grad_obj", &qpscaling_lb_norm_inf_grad_obj);

    qpscaling_scale_objective_type qpscaling_scale_objective = {{ solver_options.qpscaling_scale_objective }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qpscaling_scale_objective", &qpscaling_scale_objective);

    ocp_nlp_qpscaling_constraint_type qpscaling_scale_constraints = {{ solver_options.qpscaling_scale_constraints }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qpscaling_scale_constraints", &qpscaling_scale_constraints);

    // NLP QP tol strategy
    ocp_nlp_qp_tol_strategy_t nlp_qp_tol_strategy = {{ solver_options.nlp_qp_tol_strategy }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "nlp_qp_tol_strategy", &nlp_qp_tol_strategy);

    double nlp_qp_tol_reduction_factor = {{ solver_options.nlp_qp_tol_reduction_factor }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "nlp_qp_tol_reduction_factor", &nlp_qp_tol_reduction_factor);

    double nlp_qp_tol_safety_factor = {{ solver_options.nlp_qp_tol_safety_factor }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "nlp_qp_tol_safety_factor", &nlp_qp_tol_safety_factor);

    double nlp_qp_tol_min_stat = {{ solver_options.nlp_qp_tol_min_stat }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "nlp_qp_tol_min_stat", &nlp_qp_tol_min_stat);

    double nlp_qp_tol_min_eq = {{ solver_options.nlp_qp_tol_min_eq }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "nlp_qp_tol_min_eq", &nlp_qp_tol_min_eq);

    double nlp_qp_tol_min_ineq = {{ solver_options.nlp_qp_tol_min_ineq }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "nlp_qp_tol_min_ineq", &nlp_qp_tol_min_ineq);

    double nlp_qp_tol_min_comp = {{ solver_options.nlp_qp_tol_min_comp }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "nlp_qp_tol_min_comp", &nlp_qp_tol_min_comp);

{%- if solver_options.nlp_solver_type == "SQP" and solver_options.timeout_max_time > 0 %}
    double timeout_max_time = {{ solver_options.timeout_max_time }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "timeout_max_time", &timeout_max_time);

    ocp_nlp_timeout_heuristic_t timeout_heuristic = {{ solver_options.timeout_heuristic }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "timeout_heuristic", &timeout_heuristic);
{%- endif %}

{%- elif solver_options.nlp_solver_type == "SQP_RTI" %}
    int as_rti_iter = {{ solver_options.as_rti_iter }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "as_rti_iter", &as_rti_iter);

    int as_rti_level = {{ solver_options.as_rti_level }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "as_rti_level", &as_rti_level);

    int rti_log_residuals = {{ solver_options.rti_log_residuals }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "rti_log_residuals", &rti_log_residuals);

    int rti_log_only_available_residuals = {{ solver_options.rti_log_only_available_residuals }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "rti_log_only_available_residuals", &rti_log_only_available_residuals);
{%- endif %}

    bool with_anderson_acceleration = {{ solver_options.with_anderson_acceleration }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "with_anderson_acceleration", &with_anderson_acceleration);

    int qp_solver_iter_max = {{ solver_options.qp_solver_iter_max }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_iter_max", &qp_solver_iter_max);

{# NOTE: qp_solver tolerances must be set after NLP ones, since the setter for NLP tolerances sets the QP tolerances to the same values. #}
    {%- if solver_options.qp_solver_tol_stat %}
    double qp_solver_tol_stat = {{ solver_options.qp_solver_tol_stat }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_tol_stat", &qp_solver_tol_stat);
    {%- endif %}

    {%- if solver_options.qp_solver_tol_eq %}
    double qp_solver_tol_eq = {{ solver_options.qp_solver_tol_eq }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_tol_eq", &qp_solver_tol_eq);
    {%- endif %}

    {%- if solver_options.qp_solver_tol_ineq %}
    double qp_solver_tol_ineq = {{ solver_options.qp_solver_tol_ineq }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_tol_ineq", &qp_solver_tol_ineq);
    {%- endif %}

    {%- if solver_options.qp_solver_tol_comp %}
    double qp_solver_tol_comp = {{ solver_options.qp_solver_tol_comp }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_tol_comp", &qp_solver_tol_comp);
    {%- endif %}

    {%- if solver_options.qp_solver_warm_start %}
    int qp_solver_warm_start = {{ solver_options.qp_solver_warm_start }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_warm_start", &qp_solver_warm_start);
    {%- endif %}

    int print_level = {{ solver_options.print_level }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "print_level", &print_level);

{%- if solver_options.qp_solver is containing('PARTIAL_CONDENSING') %}
    int qp_solver_cond_ric_alg = {{ solver_options.qp_solver_cond_ric_alg }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_cond_ric_alg", &qp_solver_cond_ric_alg);
{% endif %}

{%- if solver_options.qp_solver == 'PARTIAL_CONDENSING_HPIPM' %}
    int qp_solver_ric_alg = {{ solver_options.qp_solver_ric_alg }};
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_ric_alg", &qp_solver_ric_alg);
{% endif %}


    /* Stage varying options */
    int ext_cost_num_hess;
    bool output_z_val = true;
    bool sens_algebraic_val = true;
    sim_collocation_type collocation_type;

    // set up sim_method_num_stages
    int* sim_method_num_stages = malloc(N*sizeof(int));
    {%- for j in range(end=N_horizon) %}
    sim_method_num_stages[{{ j }}] = {{ solver_options.sim_method_num_stages[j] }};
    {%- endfor %}

    // set up sim_method_num_steps
    int* sim_method_num_steps = malloc(N*sizeof(int));
    {%- for j in range(end=N_horizon) %}
    sim_method_num_steps[{{ j }}] = {{ solver_options.sim_method_num_steps[j] }};
    {%- endfor %}

    // set up sim_method_jac_reuse
    bool* sim_method_jac_reuse = malloc(N*sizeof(bool));
    {%- for j in range(end=N_horizon) %}
    sim_method_jac_reuse[{{ j }}] = (bool){{ solver_options.sim_method_jac_reuse[j] }};
    {%- endfor %}

    {%- for jj in range(end=n_phases) %}{# phases loop !#}

{%- if phases_dims[jj].nz > 0 %}
    // TODO: these options are lower level -> should be encapsulated! maybe through hessian approx option.
    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_output_z", &output_z_val);
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_sens_algebraic", &sens_algebraic_val);
    }
{%- endif %}

{%- if mocp_opts.integrator_type[jj] != "DISCRETE" %}

    // set collocation type (relevant for implicit integrators)
    collocation_type = {{ mocp_opts.collocation_type[jj] }};
    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_collocation_type", &collocation_type);

    // sim_method_newton_iter
    newton_iter_val = {{ solver_options.sim_method_newton_iter }};
    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_newton_iter", &newton_iter_val);

    // sim_method_newton_tol
    newton_tol_val = {{ solver_options.sim_method_newton_tol }};
    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_newton_tol", &newton_tol_val);

    // possibly varying: sim_method_num_steps, sim_method_num_stages, sim_method_jac_reuse
    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_num_steps", &sim_method_num_steps[i]);

    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_num_stages", &sim_method_num_stages[i]);

    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_jac_reuse", &sim_method_jac_reuse[i]);

{%- if mocp_opts.cost_discretization[jj] == "INTEGRATOR" %}
    tmp_bool = true;
    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_solver_opts_set_at_stage(nlp_config, capsule->nlp_opts, i, "dynamics_cost_computation", &tmp_bool);
        ocp_nlp_solver_opts_set_at_stage(nlp_config, capsule->nlp_opts, i, "dynamics_cost_type", &capsule->nlp_solver_plan->nlp_cost[i]);
        ocp_nlp_solver_opts_set_at_stage(nlp_config, capsule->nlp_opts, i, "cost_integrator_cost", &tmp_bool);
    }
{%- endif %}
{%- endif %}{# mocp_opts.integrator_type[jj] != "DISCRETE" #}

{%- if cost[jj].cost_type == "EXTERNAL" %}
    ext_cost_num_hess = {{ solver_options.ext_cost_num_hess }};
    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "cost_numerical_hessian", &ext_cost_num_hess);
    }
{%- endif %}
    {%- endfor %}{# for jj in range(end=n_phases) #}

{%- if cost_e.cost_type_e == "EXTERNAL" %}
    ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, N, "cost_numerical_hessian", &ext_cost_num_hess);
{%- endif %}

    // free arrays
    free(sim_method_num_steps);
    free(sim_method_num_stages);
    free(sim_method_jac_reuse);
}

/**
 * Internal function for {{ name }}_acados_create: step 7
 */
void {{ name }}_acados_create_set_nlp_out({{ name }}_solver_capsule* capsule)
{
    const int N = capsule->nlp_solver_plan->N;
    ocp_nlp_config* nlp_config = capsule->nlp_config;
    ocp_nlp_dims* nlp_dims = capsule->nlp_dims;
    ocp_nlp_out* nlp_out = capsule->nlp_out;
    ocp_nlp_in* nlp_in = capsule->nlp_in;


    int nx_max = {{ nx_max }};
    int nu_max = {{ nu_max }};

    // initialize primal solution
    double* xu0 = calloc(nx_max+nu_max, sizeof(double));
    double* x0 = xu0;
{% if dims_0.nbx_0 == dims_0.nx %}
    // initialize with x0
    {%- for item in constraints_0.lbx_0 %}
        {%- if item != 0 %}
    x0[{{ loop.index0 }}] = {{ item }};
        {%- endif %}
    {%- endfor %}
{% else %}
    // initialize with zeros
{%- endif %}

    double* u0 = xu0 + nx_max;

    for (int i = 0; i < N; i++)
    {
        // x0
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, i, "x", x0);
        // u0
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, i, "u", u0);
    }
    ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, N, "x", x0);
    free(xu0);
}



/**
 * Internal function for {{ name }}_acados_create: step 9
 */
int {{ name }}_acados_create_precompute({{ name }}_solver_capsule* capsule) {
    int status = ocp_nlp_precompute(capsule->nlp_solver, capsule->nlp_in, capsule->nlp_out);

    if (status != ACADOS_SUCCESS) {
        printf("\nocp_nlp_precompute failed!\n\n");
        exit(1);
    }

    return status;
}



int {{ name }}_acados_create_with_discretization({{ name }}_solver_capsule* capsule, int N, double* new_time_steps)
{
    // If N does not match the number of shooting intervals used for code generation, new_time_steps must be given.
    if (new_time_steps) {
        fprintf(stderr, "{{ name }}_acados_create_with_discretization: new_time_steps should be NULL " \
            "for multi-phase solver!\n");
        return 1;
    }

    // 1) create and set nlp_solver_plan; create nlp_config
    capsule->nlp_solver_plan = ocp_nlp_plan_create(N);
    {{ name }}_acados_create_set_plan(capsule->nlp_solver_plan, N);
    capsule->nlp_config = ocp_nlp_config_create(*capsule->nlp_solver_plan);

    // 2) create and set dimensions
    capsule->nlp_dims = {{ name }}_acados_create_setup_dimensions(capsule);

    // 3) create and set nlp_opts
    capsule->nlp_opts = ocp_nlp_solver_opts_create(capsule->nlp_config, capsule->nlp_dims);
    {{ name }}_acados_create_set_opts(capsule);

    // 4) create and set nlp_out
    // 4.1) nlp_out
    capsule->nlp_out = ocp_nlp_out_create(capsule->nlp_config, capsule->nlp_dims);
    // 4.2) sens_out
    capsule->sens_out = ocp_nlp_out_create(capsule->nlp_config, capsule->nlp_dims);
    {{ name }}_acados_create_set_nlp_out(capsule);

    // 5) create nlp_in
    capsule->nlp_in = ocp_nlp_in_create(capsule->nlp_config, capsule->nlp_dims);

    // 6) set default parameters in functions
    {{ name }}_acados_create_setup_functions(capsule);
    {{ name }}_acados_create_setup_nlp_in(capsule, N);
    {{ name }}_acados_create_set_default_parameters(capsule);

    // 7) create solver
    capsule->nlp_solver = ocp_nlp_solver_create(capsule->nlp_config, capsule->nlp_dims, capsule->nlp_opts, capsule->nlp_in);


    // 8) do precomputations
    int status = {{ name }}_acados_create_precompute(capsule);

    {%- if custom_update_filename != "" %}
    // Initialize custom update function
    custom_update_init_function(capsule);
    {%- endif %}

    return status;
}




int {{ name }}_acados_reset({{ name }}_solver_capsule* capsule, int reset_qp_solver_mem)
{
    // set initialization to all zeros
{# TODO: use guess values / initial state value from json instead?! #}
    const int N = capsule->nlp_solver_plan->N;
    ocp_nlp_config* nlp_config = capsule->nlp_config;
    ocp_nlp_dims* nlp_dims = capsule->nlp_dims;
    ocp_nlp_out* nlp_out = capsule->nlp_out;
    ocp_nlp_in* nlp_in = capsule->nlp_in;
    ocp_nlp_solver* nlp_solver = capsule->nlp_solver;

{%- set_global buffer_size = 0 %}
{%- for jj in range(end=n_phases) %}{# phases loop !#}
    {%- set buffer_sizes = [
            phases_dims[jj].nx,
            phases_dims[jj].nu,
            phases_dims[jj].nz,
            2 * phases_dims[jj].ns,
            2 * phases_dims[jj].ns_0,
            2 * phases_dims[jj].ns_e,
            phases_dims[jj].nbx_0,
            phases_dims[jj].nbx,
            phases_dims[jj].nbx_e,
            phases_dims[jj].nbu,
            phases_dims[jj].ng,
            phases_dims[jj].ng_e,
            phases_dims[jj].nh,
            phases_dims[jj].nh_0,
            phases_dims[jj].nh_e,
            phases_dims[jj].nphi,
            phases_dims[jj].nphi_0,
            phases_dims[jj].nphi_e,
        ]
    -%}
{%- for b in buffer_sizes %}
{%- set_global buffer_size = buffer_size + b %}
{%- endfor %}
{%- endfor %}

    double* buffer = calloc({{ buffer_size }}, sizeof(double));


{%- for jj in range(end=n_phases) %}{# phases loop !#}
    // Reset stage {{ jj }}
    for (int i = {{ start_idx[jj] }}; i < {{ end_idx[jj] }}; i++)
    {
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, i, "x", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, i, "u", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, i, "sl", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, i, "su", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, i, "lam", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, i, "z", buffer);
        if (i<N)
        {
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, i, "pi", buffer);
        {%- if mocp_opts.integrator_type[jj] == "IRK" %}
            ocp_nlp_set(nlp_solver, i, "xdot_guess", buffer);
            ocp_nlp_set(nlp_solver, i, "z_guess", buffer);
        {%- elif mocp_opts.integrator_type[jj] == "LIFTED_IRK" %}
            ocp_nlp_set(nlp_solver, i, "xdot_guess", buffer);
        {%- elif mocp_opts.integrator_type[jj] == "GNSF" %}
            ocp_nlp_set(nlp_solver, i, "gnsf_phi_guess", buffer);
        {%- endif %}
        }
    }
{%- endfor %}

{%- if solver_options.qp_solver == 'PARTIAL_CONDENSING_HPIPM' %}
    // get qp_status: if NaN -> reset memory
    int qp_status;
    ocp_nlp_get(capsule->nlp_solver, "qp_status", &qp_status);
    if (reset_qp_solver_mem || (qp_status == 3))
    {
        // printf("\nin reset qp_status %d -> resetting QP memory\n", qp_status);
        ocp_nlp_solver_reset_qp_memory(nlp_solver, nlp_in, nlp_out);
    }
{%- endif %}

    free(buffer);
    return 0;
}




int {{ name }}_acados_update_params({{ name }}_solver_capsule* capsule, int stage, double *p, int np)
{
    int solver_status = 0;

    {%- for jj in range(end=n_phases) %}{# phases loop !#}
    if (stage >= {{ start_idx[jj] }} && stage < {{ end_idx[jj] }})
    {
        if (NP_{{ jj }} != np)
        {
            printf("acados_update_params: trying to set %i parameters at stage %i."
                " Parameters should be of length %i. Exiting.\n", np, stage, NP_{{ jj }});
            exit(1);
        }
    }
    {% endfor %}
    ocp_nlp_in_set(capsule->nlp_config, capsule->nlp_dims, capsule->nlp_in, stage, "parameter_values", p);

    return solver_status;
}


int {{ name }}_acados_update_params_sparse({{ name }}_solver_capsule * capsule, int stage, int *idx, double *p, int n_update)
{
    ocp_nlp_in_set_params_sparse(capsule->nlp_config, capsule->nlp_dims, capsule->nlp_in, stage, idx, p, n_update);

    return 0;
}


int {{ name }}_acados_set_p_global_and_precompute_dependencies({{ name }}_solver_capsule* capsule, double* data, int data_len)
{
{% if dims_0.n_global_data > 0 %}
    external_function_casadi* fun = &capsule->p_global_precompute_fun;
    fun->args[0] = data;
    int np_global = {{ dims_0.np_global }};

    if (data_len != np_global)
    {
        printf("{{ name }}_acados_set_p_global_and_precompute_dependencies: np_global = %d should match data_len = %d. Exiting.\n", np_global, data_len);
        exit(1);
    }

    ocp_nlp_in *in = {{ name }}_acados_get_nlp_in(capsule);
    fun->res[0] = in->global_data;

    fun->casadi_fun((const double **) fun->args, fun->res, fun->int_work, fun->float_work, NULL);

{%- else %}
    // printf("No global_data, {{ name }}_acados_set_p_global_and_precompute_dependencies does nothing.\n");
{%- endif %}
    return 0;
}



int {{ name }}_acados_solve({{ name }}_solver_capsule* capsule)
{
    // solve NLP
    int solver_status = ocp_nlp_solve(capsule->nlp_solver, capsule->nlp_in, capsule->nlp_out);

    return solver_status;
}


int {{ name }}_acados_setup_qp_matrices_and_factorize({{ name }}_solver_capsule* capsule)
{
    // solve NLP
    int solver_status = ocp_nlp_setup_qp_matrices_and_factorize(capsule->nlp_solver, capsule->nlp_in, capsule->nlp_out);

    return solver_status;
}



void {{ name }}_acados_print_stats({{ name }}_solver_capsule* capsule)
{
    int sqp_iter, stat_m, stat_n, tmp_int;
    ocp_nlp_get(capsule->nlp_solver, "sqp_iter", &sqp_iter);
    ocp_nlp_get(capsule->nlp_solver, "stat_n", &stat_n);
    ocp_nlp_get(capsule->nlp_solver, "stat_m", &stat_m);

{% set stat_n_max = 12 %}
    double stat[{{ solver_options.nlp_solver_max_iter * stat_n_max }}];
    ocp_nlp_get(capsule->nlp_solver, "statistics", stat);

    int nrow = sqp_iter+1 < stat_m ? sqp_iter+1 : stat_m;

    printf("iter\tres_stat\tres_eq\t\tres_ineq\tres_comp\tqp_stat\tqp_iter\talpha");
    if (stat_n > 8)
        printf("\t\tqp_res_stat\tqp_res_eq\tqp_res_ineq\tqp_res_comp");
    printf("\n");

{%- if solver_options.nlp_solver_type == "SQP" %}

    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < stat_n + 1; j++)
        {
            if (j == 0 || j == 5 || j == 6)
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
{% else %}
    printf("iter\tqp_stat\tqp_iter\n");
    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < stat_n + 1; j++)
        {
            tmp_int = (int) stat[i + j * nrow];
            printf("%d\t", tmp_int);
        }
        printf("\n");
    }
{%- endif %}
}




int {{ name }}_acados_free({{ name }}_solver_capsule* capsule)
{
    // before destroying, keep some info
    const int N = capsule->nlp_solver_plan->N;
    {%- if custom_update_filename != "" %}
    custom_update_terminate_function(capsule);
    {%- endif %}
    // free memory
    ocp_nlp_solver_opts_destroy(capsule->nlp_opts);
    ocp_nlp_in_destroy(capsule->nlp_in);
    ocp_nlp_out_destroy(capsule->nlp_out);
    ocp_nlp_out_destroy(capsule->sens_out);
    ocp_nlp_solver_destroy(capsule->nlp_solver);
    ocp_nlp_dims_destroy(capsule->nlp_dims);
    ocp_nlp_config_destroy(capsule->nlp_config);
    ocp_nlp_plan_destroy(capsule->nlp_solver_plan);

    /* free external function */
    // initial node
{%- if cost[0].cost_type_0 == "NONLINEAR_LS" %}
    external_function_external_param_casadi_free(&capsule->cost_y_0_fun);
    external_function_external_param_casadi_free(&capsule->cost_y_0_fun_jac_ut_xt);
    {%- if solver_options.hessian_approx == "EXACT" %}
    external_function_external_param_casadi_free(&capsule->cost_y_0_hess);
    {%- endif %}
{%- elif cost[0].cost_type_0 == "CONVEX_OVER_NONLINEAR" %}
    external_function_external_param_casadi_free(&capsule->conl_cost_0_fun);
    external_function_external_param_casadi_free(&capsule->conl_cost_0_fun_jac_hess);
{%- elif cost[0].cost_type_0 == "EXTERNAL" %}
    external_function_external_param_{{ cost[0].cost_ext_fun_type_0 }}_free(&capsule->ext_cost_0_fun);
    external_function_external_param_{{ cost[0].cost_ext_fun_type_0 }}_free(&capsule->ext_cost_0_fun_jac);
    external_function_external_param_{{ cost[0].cost_ext_fun_type_0 }}_free(&capsule->ext_cost_0_fun_jac_hess);
{%- endif %}

{%- if constraints_0.constr_type_0 == "BGH" and dims_0.nh_0 > 0 %}
    external_function_external_param_casadi_free(&capsule->nl_constr_h_0_fun_jac);
    external_function_external_param_casadi_free(&capsule->nl_constr_h_0_fun);
{%- if solver_options.hessian_approx == "EXACT" %}
    external_function_external_param_casadi_free(&capsule->nl_constr_h_0_fun_jac_hess);
{%- endif %}
{%- if solver_options.with_solution_sens_wrt_params %}
    external_function_external_param_casadi_free(&capsule->nl_constr_h_0_jac_p_hess_xu_p);
{%- endif %}
{%- if solver_options.with_value_sens_wrt_params %}
    external_function_external_param_casadi_free(&capsule->nl_constr_h_0_adj_p);
{%- endif %}
{%- elif constraints_0.constr_type_0 == "BGP" and dims_0.nphi_0 > 0 %}
    external_function_external_param_casadi_free(&capsule->phi_0_constraint_fun);
    external_function_external_param_casadi_free(&capsule->phi_0_constraint_fun_jac_hess);
{%- endif %}




{%- for jj in range(end=n_phases) %}{# phases loop !#}
    /* Path phase {jj} */
    // dynamics
{%- if mocp_opts.integrator_type[jj] == "IRK" %}
    for (int i_fun = 0; i_fun < {{ end_idx[jj] - start_idx[jj] }}; i_fun++)
    {
        external_function_external_param_{{ model[jj].dyn_ext_fun_type }}_free(&capsule->impl_dae_fun_{{ jj }}[i_fun]);
        external_function_external_param_{{ model[jj].dyn_ext_fun_type }}_free(&capsule->impl_dae_fun_jac_x_xdot_z_{{ jj }}[i_fun]);
        external_function_external_param_{{ model[jj].dyn_ext_fun_type }}_free(&capsule->impl_dae_jac_x_xdot_u_z_{{ jj }}[i_fun]);
    {%- if solver_options.hessian_approx == "EXACT" %}
        external_function_external_param_{{ model[jj].dyn_ext_fun_type }}_free(&capsule->impl_dae_hess_{{ jj }}[i_fun]);
    {%- endif %}
    }
    free(capsule->impl_dae_fun_{{ jj }});
    free(capsule->impl_dae_fun_jac_x_xdot_z_{{ jj }});
    free(capsule->impl_dae_jac_x_xdot_u_z_{{ jj }});
    {%- if solver_options.hessian_approx == "EXACT" %}
    free(capsule->impl_dae_hess_{{ jj }});
    {%- endif %}

{%- elif mocp_opts.integrator_type[jj] == "LIFTED_IRK" %}
    for (int i_fun = 0; i_fun < {{ end_idx[jj] - start_idx[jj] }}; i_fun++)
    {
        external_function_external_param_{{ model[jj].dyn_ext_fun_type }}_free(&capsule->impl_dae_fun_{{ jj }}[i_fun]);
        external_function_external_param_{{ model[jj].dyn_ext_fun_type }}_free(&capsule->impl_dae_fun_jac_x_xdot_u_{{ jj }}[i_fun]);
    }
    free(capsule->impl_dae_fun_{{ jj }});
    free(capsule->impl_dae_fun_jac_x_xdot_u_{{ jj }});

{%- elif mocp_opts.integrator_type[jj] == "ERK" %}
    for (int i_fun = 0; i_fun < {{ end_idx[jj] - start_idx[jj] }}; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->expl_vde_forw_{{ jj }}[i_fun]);
        external_function_external_param_casadi_free(&capsule->expl_vde_adj_{{ jj }}[i_fun]);
        external_function_external_param_casadi_free(&capsule->expl_ode_fun_{{ jj }}[i_fun]);
    {%- if solver_options.hessian_approx == "EXACT" %}
        external_function_external_param_casadi_free(&capsule->expl_ode_hess_{{ jj }}[i_fun]);
    {%- endif %}
    }
    free(capsule->expl_vde_forw_{{ jj }});
    free(capsule->expl_vde_adj_{{ jj }});
    free(capsule->expl_ode_fun_{{ jj }});
    {%- if solver_options.hessian_approx == "EXACT" %}
    free(capsule->expl_ode_hess_{{ jj }});
    {%- endif %}

{%- elif mocp_opts.integrator_type[jj] == "GNSF" %}
    for (int i_fun = 0; i_fun < {{ end_idx[jj] - start_idx[jj] }}; i_fun++)
    {
        {% if model[jj].gnsf_purely_linear != 1 %}
        external_function_external_param_casadi_free(&capsule->gnsf_phi_fun_{{ jj }}[i_fun]);
        external_function_external_param_casadi_free(&capsule->gnsf_phi_fun_jac_y_{{ jj }}[i_fun]);
        external_function_external_param_casadi_free(&capsule->gnsf_phi_jac_y_uhat_{{ jj }}[i_fun]);
        {% if model[jj].gnsf_nontrivial_f_LO == 1 %}
        external_function_external_param_casadi_free(&capsule->gnsf_f_lo_jac_x1_x1dot_u_z_{{ jj }}[i_fun]);
        {%- endif %}
        {%- endif %}
        external_function_external_param_casadi_free(&capsule->gnsf_get_matrices_fun_{{ jj }}[i_fun]);
    }
  {% if model[jj].gnsf_purely_linear != 1 %}
    free(capsule->gnsf_phi_fun_{{ jj }});
    free(capsule->gnsf_phi_fun_jac_y_{{ jj }});
    free(capsule->gnsf_phi_jac_y_uhat_{{ jj }});
  {% if model[jj].gnsf_nontrivial_f_LO == 1 %}
    free(capsule->gnsf_f_lo_jac_x1_x1dot_u_z_{{ jj }});
  {%- endif %}
  {%- endif %}
    free(capsule->gnsf_get_matrices_fun_{{ jj }});
{%- elif mocp_opts.integrator_type[jj] == "DISCRETE" %}
    for (int i_fun = 0; i_fun < {{ end_idx[jj] - start_idx[jj] }}; i_fun++)
    {
        external_function_external_param_{{ model[jj].dyn_ext_fun_type }}_free(&capsule->discr_dyn_phi_fun_{{ jj }}[i_fun]);
        external_function_external_param_{{ model[jj].dyn_ext_fun_type }}_free(&capsule->discr_dyn_phi_fun_jac_ut_xt_{{ jj }}[i_fun]);
    {%- if solver_options.hessian_approx == "EXACT" %}
        external_function_external_param_{{ model[jj].dyn_ext_fun_type }}_free(&capsule->discr_dyn_phi_fun_jac_ut_xt_hess_{{ jj }}[i_fun]);
    {%- endif %}
    }
    free(capsule->discr_dyn_phi_fun_{{ jj }});
    free(capsule->discr_dyn_phi_fun_jac_ut_xt_{{ jj }});
    {%- if solver_options.hessian_approx == "EXACT" %}
    free(capsule->discr_dyn_phi_fun_jac_ut_xt_hess_{{ jj }});
    {%- endif %}
{%- endif %}


{%- if cost[jj].cost_type == "NONLINEAR_LS" %}
    for (int i_fun = 0; i_fun < {{ end_idx[jj] - cost_start_idx[jj] }}; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->cost_y_fun_{{ jj }}[i_fun]);
        external_function_external_param_casadi_free(&capsule->cost_y_fun_jac_ut_xt_{{ jj }}[i_fun]);
        {%- if solver_options.hessian_approx == "EXACT" %}
        external_function_external_param_casadi_free(&capsule->cost_y_hess_{{ jj }}[i_fun]);
        {%- endif %}
    }
    free(capsule->cost_y_fun_{{ jj }});
    free(capsule->cost_y_fun_jac_ut_xt_{{ jj }});
    {%- if solver_options.hessian_approx == "EXACT" %}
    free(capsule->cost_y_hess_{{ jj }});
    {%- endif %}

{%- elif cost[jj].cost_type == "CONVEX_OVER_NONLINEAR" %}
    for (int i_fun = 0; i_fun < {{ end_idx[jj] - cost_start_idx[jj] }}; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->conl_cost_fun_{{ jj }}[i_fun]);
        external_function_external_param_casadi_free(&capsule->conl_cost_fun_jac_hess_{{ jj }}[i_fun]);
    }
    free(capsule->conl_cost_fun_{{ jj }});
    free(capsule->conl_cost_fun_jac_hess_{{ jj }});
{%- elif cost[jj].cost_type == "EXTERNAL" %}
    for (int i_fun = 0; i_fun < {{ end_idx[jj] - cost_start_idx[jj] }}; i_fun++)
    {
        external_function_external_param_{{ cost[jj].cost_ext_fun_type }}_free(&capsule->ext_cost_fun_{{ jj }}[i_fun]);
        external_function_external_param_{{ cost[jj].cost_ext_fun_type }}_free(&capsule->ext_cost_fun_jac_{{ jj }}[i_fun]);
        external_function_external_param_{{ cost[jj].cost_ext_fun_type }}_free(&capsule->ext_cost_fun_jac_hess_{{ jj }}[i_fun]);
    }
    free(capsule->ext_cost_fun_{{ jj }});
    free(capsule->ext_cost_fun_jac_{{ jj }});
    free(capsule->ext_cost_fun_jac_hess_{{ jj }});
{%- endif %}

    // constraints
{%- if constraints[jj].constr_type == "BGH" and phases_dims[jj].nh > 0 %}
    for (int i_fun = 0; i_fun < {{ end_idx[jj] - cost_start_idx[jj] }}; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->nl_constr_h_fun_jac_{{ jj }}[i_fun]);
        external_function_external_param_casadi_free(&capsule->nl_constr_h_fun_{{ jj }}[i_fun]);
  {%- if solver_options.hessian_approx == "EXACT" %}
        external_function_external_param_casadi_free(&capsule->nl_constr_h_fun_jac_hess_{{ jj }}[i_fun]);
  {%- endif %}
  {%- if solver_options.with_solution_sens_wrt_params %}
        external_function_external_param_casadi_free(&capsule->nl_constr_h_jac_p_hess_xu_p_{{ jj }}[i_fun]);
  {%- endif %}
  {%- if solver_options.with_solution_sens_wrt_params %}
        external_function_external_param_casadi_free(&capsule->nl_constr_h_adj_p_{{ jj }}[i_fun]);
  {%- endif %}
    }
    free(capsule->nl_constr_h_fun_jac_{{ jj }});
    free(capsule->nl_constr_h_fun_{{ jj }});
  {%- if solver_options.hessian_approx == "EXACT" %}
    free(capsule->nl_constr_h_fun_jac_hess_{{ jj }});
  {%- endif %}
  {%- if solver_options.with_solution_sens_wrt_params %}
    free(capsule->nl_constr_h_jac_p_hess_xu_p_{{ jj }});
  {%- endif %}
  {%- if solver_options.with_solution_sens_wrt_params %}
    free(capsule->nl_constr_h_adj_p_{{ jj }});
  {%- endif %}

{%- elif constraints[jj].constr_type == "BGP" and phases_dims[jj].nphi > 0 %}
    for (int i_fun = 0; i_fun < {{ end_idx[jj] - cost_start_idx[jj] }}; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->phi_constraint_fun_{{ jj }}[i_fun]);
        external_function_external_param_casadi_free(&capsule->phi_constraint_fun_jac_hess_{{ jj }}[i_fun]);
    }
    free(capsule->phi_constraint_fun_{{ jj }});
    free(capsule->phi_constraint_fun_jac_hess_{{ jj }});
{%- endif %}
    {%- endfor %}{# for jj in range(end=n_phases) #}


    /* Terminal node */
{%- if constraints_e.constr_type_e == "BGH" and dims_e.nh_e > 0 %}
    external_function_external_param_casadi_free(&capsule->nl_constr_h_e_fun_jac);
    external_function_external_param_casadi_free(&capsule->nl_constr_h_e_fun);
{%- if solver_options.hessian_approx == "EXACT" %}
    external_function_external_param_casadi_free(&capsule->nl_constr_h_e_fun_jac_hess);
{%- endif %}
{%- if solver_options.with_solution_sens_wrt_params %}
    external_function_external_param_casadi_free(&capsule->nl_constr_h_e_jac_p_hess_xu_p);
{%- endif %}
{%- if solver_options.with_value_sens_wrt_params %}
    external_function_external_param_casadi_free(&capsule->nl_constr_h_e_adj_p);
{%- endif %}
{%- elif constraints_e.constr_type_e == "BGP" and dims_e.nphi_e > 0 %}
    external_function_external_param_casadi_free(&capsule->phi_e_constraint_fun);
    external_function_external_param_casadi_free(&capsule->phi_e_constraint_fun_jac_hess);
{%- endif %}

{%- if cost_e.cost_type_e == "NONLINEAR_LS" %}
    external_function_external_param_casadi_free(&capsule->cost_y_e_fun);
    external_function_external_param_casadi_free(&capsule->cost_y_e_fun_jac_ut_xt);
    {%- if solver_options.hessian_approx == "EXACT" %}
    external_function_external_param_casadi_free(&capsule->cost_y_e_hess);
    {%- endif %}
{%- elif cost_e.cost_type_e == "CONVEX_OVER_NONLINEAR" %}
    external_function_external_param_casadi_free(&capsule->conl_cost_e_fun);
    external_function_external_param_casadi_free(&capsule->conl_cost_e_fun_jac_hess);
{%- elif cost_e.cost_type_e == "EXTERNAL" %}
    external_function_external_param_{{ cost_e.cost_ext_fun_type_e }}_free(&capsule->ext_cost_e_fun);
    external_function_external_param_{{ cost_e.cost_ext_fun_type_e }}_free(&capsule->ext_cost_e_fun_jac);
    external_function_external_param_{{ cost_e.cost_ext_fun_type_e }}_free(&capsule->ext_cost_e_fun_jac_hess);
{%- endif %}

{% if phases_dims[0].n_global_data > 0 %}
    external_function_casadi_free(&capsule->p_global_precompute_fun);
{%- endif %}

    return 0;
}




int {{ name }}_acados_custom_update({{ name }}_solver_capsule* capsule, double* data, int data_len)
{
{%- if custom_update_filename == "" %}
    (void)capsule;
    (void)data;
    (void)data_len;
    printf("\ndummy function that can be called in between solver calls to update parameters or numerical data efficiently in C.\n");
    printf("nothing set yet..\n");
    return 1;
{% else %}
    custom_update_function(capsule, data, data_len);
{%- endif %}
}


ocp_nlp_in *{{ name }}_acados_get_nlp_in({{ name }}_solver_capsule* capsule) { return capsule->nlp_in; }
ocp_nlp_out *{{ name }}_acados_get_nlp_out({{ name }}_solver_capsule* capsule) { return capsule->nlp_out; }
ocp_nlp_out *{{ name }}_acados_get_sens_out({{ name }}_solver_capsule* capsule) { return capsule->sens_out; }
ocp_nlp_solver *{{ name }}_acados_get_nlp_solver({{ name }}_solver_capsule* capsule) { return capsule->nlp_solver; }
ocp_nlp_config *{{ name }}_acados_get_nlp_config({{ name }}_solver_capsule* capsule) { return capsule->nlp_config; }
void *{{ name }}_acados_get_nlp_opts({{ name }}_solver_capsule* capsule) { return capsule->nlp_opts; }
ocp_nlp_dims *{{ name }}_acados_get_nlp_dims({{ name }}_solver_capsule* capsule) { return capsule->nlp_dims; }
ocp_nlp_plan_t *{{ name }}_acados_get_nlp_plan({{ name }}_solver_capsule* capsule) { return capsule->nlp_solver_plan; }
