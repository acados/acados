/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */

#ifndef ACADOS_SIM_{{ name }}_H_
#define ACADOS_SIM_{{ name }}_H_

#include "acados_c/sim_interface.h"
#include "acados_c/external_function_interface.h"

#define {{ name | upper }}_NX     {{ dims.nx }}
#define {{ name | upper }}_NZ     {{ dims.nz }}
#define {{ name | upper }}_NU     {{ dims.nu }}
#define {{ name | upper }}_NP     {{ dims.np }}

#ifdef __cplusplus
extern "C" {
#endif


// ** capsule for solver data **
typedef struct {{ name }}_sim_solver_capsule
{
    // acados objects
    sim_in *acados_sim_in;
    sim_out *acados_sim_out;
    sim_solver *acados_sim_solver;
    sim_opts *acados_sim_opts;
    sim_config *acados_sim_config;
    void *acados_sim_dims;
    void *acados_sim_mem;

    /* external functions */
    // ERK
    external_function_param_{{ model.dyn_ext_fun_type }} * sim_expl_vde_forw;
    external_function_param_{{ model.dyn_ext_fun_type }} * sim_vde_adj_casadi;
    external_function_param_{{ model.dyn_ext_fun_type }} * sim_expl_ode_fun_casadi;
    external_function_param_{{ model.dyn_ext_fun_type }} * sim_expl_ode_hess;
    external_function_param_{{ model.dyn_ext_fun_type }} * sim_expl_vde_forw_p;

    // IRK
    external_function_param_{{ model.dyn_ext_fun_type }} * sim_impl_dae_fun;
    external_function_param_{{ model.dyn_ext_fun_type }} * sim_impl_dae_fun_jac_x_xdot_z;
    external_function_param_{{ model.dyn_ext_fun_type }} * sim_impl_dae_jac_x_xdot_u_z;
    external_function_param_{{ model.dyn_ext_fun_type }} * sim_impl_dae_hess;
    external_function_param_{{ model.dyn_ext_fun_type }} * sim_impl_dae_jac_p;

    // GNSF
    external_function_param_{{ model.dyn_ext_fun_type }} * sim_gnsf_phi_fun;
    external_function_param_{{ model.dyn_ext_fun_type }} * sim_gnsf_phi_fun_jac_y;
    external_function_param_{{ model.dyn_ext_fun_type }} * sim_gnsf_phi_jac_y_uhat;
    external_function_param_{{ model.dyn_ext_fun_type }} * sim_gnsf_f_lo_jac_x1_x1dot_u_z;
    external_function_param_{{ model.dyn_ext_fun_type }} * sim_gnsf_get_matrices_fun;

} {{ name }}_sim_solver_capsule;


ACADOS_SYMBOL_EXPORT int {{ name }}_acados_sim_create({{ name }}_sim_solver_capsule *capsule);
ACADOS_SYMBOL_EXPORT int {{ name }}_acados_sim_solve({{ name }}_sim_solver_capsule *capsule);
{% if solver_options.with_batch_functionality %}
ACADOS_SYMBOL_EXPORT void {{ name }}_acados_sim_batch_solve({{ name }}_sim_solver_capsule **capsules, int N_batch, int num_threads_in_batch_solve);
{% endif %}
ACADOS_SYMBOL_EXPORT int {{ name }}_acados_sim_free({{ name }}_sim_solver_capsule *capsule);
ACADOS_SYMBOL_EXPORT int {{ name }}_acados_sim_update_params({{ name }}_sim_solver_capsule *capsule, double *value, int np);

ACADOS_SYMBOL_EXPORT sim_config * {{ name }}_acados_get_sim_config({{ name }}_sim_solver_capsule *capsule);
ACADOS_SYMBOL_EXPORT sim_in * {{ name }}_acados_get_sim_in({{ name }}_sim_solver_capsule *capsule);
ACADOS_SYMBOL_EXPORT sim_out * {{ name }}_acados_get_sim_out({{ name }}_sim_solver_capsule *capsule);
ACADOS_SYMBOL_EXPORT void * {{ name }}_acados_get_sim_dims({{ name }}_sim_solver_capsule *capsule);
ACADOS_SYMBOL_EXPORT sim_opts * {{ name }}_acados_get_sim_opts({{ name }}_sim_solver_capsule *capsule);
ACADOS_SYMBOL_EXPORT sim_solver * {{ name }}_acados_get_sim_solver({{ name }}_sim_solver_capsule *capsule);
ACADOS_SYMBOL_EXPORT void * {{ name }}_acados_get_sim_mem({{ name }}_sim_solver_capsule *capsule);

ACADOS_SYMBOL_EXPORT {{ name }}_sim_solver_capsule * {{ name }}_acados_sim_solver_create_capsule(void);
ACADOS_SYMBOL_EXPORT int {{ name }}_acados_sim_solver_free_capsule({{ name }}_sim_solver_capsule *capsule);

#ifdef __cplusplus
}
#endif

#endif  // ACADOS_SIM_{{ name }}_H_
