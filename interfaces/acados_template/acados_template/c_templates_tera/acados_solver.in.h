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

#ifndef ACADOS_SOLVER_{{ model.name }}_H_
#define ACADOS_SOLVER_{{ model.name }}_H_

#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

#define {{ model.name | upper }}_NX     {{ dims.nx }}
#define {{ model.name | upper }}_NZ     {{ dims.nz }}
#define {{ model.name | upper }}_NU     {{ dims.nu }}
#define {{ model.name | upper }}_NP     {{ dims.np }}
#define {{ model.name | upper }}_NBX    {{ dims.nbx }}
#define {{ model.name | upper }}_NBX0   {{ dims.nbx_0 }}
#define {{ model.name | upper }}_NBU    {{ dims.nbu }}
#define {{ model.name | upper }}_NSBX   {{ dims.nsbx }}
#define {{ model.name | upper }}_NSBU   {{ dims.nsbu }}
#define {{ model.name | upper }}_NSH    {{ dims.nsh }}
#define {{ model.name | upper }}_NSG    {{ dims.nsg }}
#define {{ model.name | upper }}_NSPHI  {{ dims.nsphi }}
#define {{ model.name | upper }}_NSHN   {{ dims.nsh_e }}
#define {{ model.name | upper }}_NSGN   {{ dims.nsg_e }}
#define {{ model.name | upper }}_NSPHIN {{ dims.nsphi_e }}
#define {{ model.name | upper }}_NSBXN  {{ dims.nsbx_e }}
#define {{ model.name | upper }}_NS     {{ dims.ns }}
#define {{ model.name | upper }}_NSN    {{ dims.ns_e }}
#define {{ model.name | upper }}_NG     {{ dims.ng }}
#define {{ model.name | upper }}_NBXN   {{ dims.nbx_e }}
#define {{ model.name | upper }}_NGN    {{ dims.ng_e }}
#define {{ model.name | upper }}_NY0    {{ dims.ny_0 }}
#define {{ model.name | upper }}_NY     {{ dims.ny }}
#define {{ model.name | upper }}_NYN    {{ dims.ny_e }}
#define {{ model.name | upper }}_N      {{ dims.N }}
#define {{ model.name | upper }}_NH     {{ dims.nh }}
#define {{ model.name | upper }}_NPHI   {{ dims.nphi }}
#define {{ model.name | upper }}_NHN    {{ dims.nh_e }}
#define {{ model.name | upper }}_NPHIN  {{ dims.nphi_e }}
#define {{ model.name | upper }}_NR     {{ dims.nr }}

#ifdef __cplusplus
extern "C" {
#endif

// ** capsule for solver data **
typedef struct nlp_solver_capsule
{
    // acados objects
    ocp_nlp_in *nlp_in;
    ocp_nlp_out *nlp_out;
    ocp_nlp_solver *nlp_solver;
    void *nlp_opts;
    ocp_nlp_plan *nlp_solver_plan;
    ocp_nlp_config *nlp_config;
    ocp_nlp_dims *nlp_dims;

    // number of expected runtime parameters
    unsigned int nlp_np;

    /* external functions */
    // dynamics
    external_function_param_casadi *forw_vde_casadi;
    external_function_param_casadi *expl_ode_fun;
    external_function_param_casadi *hess_vde_casadi;
    external_function_param_casadi *impl_dae_fun;
    external_function_param_casadi *impl_dae_fun_jac_x_xdot_z;
    external_function_param_casadi *impl_dae_jac_x_xdot_u_z;
    external_function_param_casadi *impl_dae_fun_jac_x_xdot_u;
    external_function_param_casadi *impl_dae_hess;
    external_function_param_casadi *gnsf_phi_fun;
    external_function_param_casadi *gnsf_phi_fun_jac_y;
    external_function_param_casadi *gnsf_phi_jac_y_uhat;
    external_function_param_casadi *gnsf_f_lo_jac_x1_x1dot_u_z;
    external_function_param_casadi *gnsf_get_matrices_fun;
    external_function_param_{{ model.dyn_ext_fun_type }} *discr_dyn_phi_fun;
    external_function_param_{{ model.dyn_ext_fun_type }} *discr_dyn_phi_fun_jac_ut_xt;
    external_function_param_{{ model.dyn_ext_fun_type }} *discr_dyn_phi_fun_jac_ut_xt_hess;

    // cost
    external_function_param_casadi *cost_y_fun;
    external_function_param_casadi *cost_y_fun_jac_ut_xt;
    external_function_param_casadi *cost_y_hess;
    external_function_param_{{ cost.cost_ext_fun_type }} *ext_cost_fun;
    external_function_param_{{ cost.cost_ext_fun_type }} *ext_cost_fun_jac;
    external_function_param_{{ cost.cost_ext_fun_type }} *ext_cost_fun_jac_hess;

    external_function_param_casadi cost_y_0_fun;
    external_function_param_casadi cost_y_0_fun_jac_ut_xt;
    external_function_param_casadi cost_y_0_hess;
    external_function_param_{{ cost.cost_ext_fun_type_0 }} ext_cost_0_fun;
    external_function_param_{{ cost.cost_ext_fun_type_0 }} ext_cost_0_fun_jac;
    external_function_param_{{ cost.cost_ext_fun_type_0 }} ext_cost_0_fun_jac_hess;

    external_function_param_casadi cost_y_e_fun;
    external_function_param_casadi cost_y_e_fun_jac_ut_xt;
    external_function_param_casadi cost_y_e_hess;
    external_function_param_{{ cost.cost_ext_fun_type_e }} ext_cost_e_fun;
    external_function_param_{{ cost.cost_ext_fun_type_e }} ext_cost_e_fun_jac;
    external_function_param_{{ cost.cost_ext_fun_type_e }} ext_cost_e_fun_jac_hess;

    // constraints
    external_function_param_casadi *phi_constraint;
    external_function_param_casadi *nl_constr_h_fun_jac;
    external_function_param_casadi *nl_constr_h_fun;
    external_function_param_casadi *nl_constr_h_fun_jac_hess;

    external_function_param_casadi phi_e_constraint;
    external_function_param_casadi nl_constr_h_e_fun_jac;
    external_function_param_casadi nl_constr_h_e_fun;
    external_function_param_casadi nl_constr_h_e_fun_jac_hess;
} nlp_solver_capsule;

nlp_solver_capsule * {{ model.name }}_acados_create_capsule(void);
int {{ model.name }}_acados_free_capsule(nlp_solver_capsule *capsule);

int {{ model.name }}_acados_create(nlp_solver_capsule * capsule);
/**
 * More generic acados_create function, since the number of stages usually is independent of the model and its
 * constraints. If the number of stages does not match the one from code-export at least a constant step_time must be
 * given to be filled in. If const_step_time=NULL, the code-exported version is used internally.
 */
int {{ model.name }}_acados_create_w_stages(nlp_solver_capsule * capsule, int n_stages, double* const_step_time);
int {{ model.name }}_acados_update_params(nlp_solver_capsule * capsule, int stage, double *value, int np);
int {{ model.name }}_acados_solve(nlp_solver_capsule * capsule);
int {{ model.name }}_acados_free(nlp_solver_capsule * capsule);
void {{ model.name }}_acados_print_stats(nlp_solver_capsule * capsule);

ocp_nlp_in *{{ model.name }}_acados_get_nlp_in(nlp_solver_capsule * capsule);
ocp_nlp_out *{{ model.name }}_acados_get_nlp_out(nlp_solver_capsule * capsule);
ocp_nlp_solver *{{ model.name }}_acados_get_nlp_solver(nlp_solver_capsule * capsule);
ocp_nlp_config *{{ model.name }}_acados_get_nlp_config(nlp_solver_capsule * capsule);
void *{{ model.name }}_acados_get_nlp_opts(nlp_solver_capsule * capsule);
ocp_nlp_dims *{{ model.name }}_acados_get_nlp_dims(nlp_solver_capsule * capsule);
ocp_nlp_plan *{{ model.name }}_acados_get_nlp_plan(nlp_solver_capsule * capsule);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SOLVER_{{ model.name }}_H_
