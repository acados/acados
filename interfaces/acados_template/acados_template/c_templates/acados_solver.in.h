#ifndef ACADOS_SOLVER_{{ra.model_name.upper()}}_H_
#define ACADOS_SOLVER_{{ra.model_name.upper()}}_H_

#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

ocp_nlp_in * nlp_in;
ocp_nlp_out * nlp_out;
ocp_nlp_solver * nlp_solver;
void * nlp_opts;
ocp_nlp_solver_plan * nlp_solver_plan;
ocp_nlp_solver_config * nlp_config;
ocp_nlp_dims * nlp_dims;
{% if ra.solver_config.integrator_type == 'ERK': %}
{% if ra.dims.np < 1: %}
external_function_casadi * forw_vde_casadi;
{% else: %}
external_function_param_casadi * forw_vde_casadi;
{% endif %}
{% if ra.solver_config.hessian_approx == 'EXACT': %} 
{% if ra.dims.np < 1: %}
external_function_casadi * hess_vde_casadi;
{% else: %}
external_function_param_casadi * hess_vde_casadi;
{% endif %}
{% endif %}
{% elif ra.solver_config.integrator_type == 'IRK': %}
{% if ra.dims.np < 1: %}
external_function_casadi * impl_dae_fun;
external_function_casadi * impl_dae_fun_jac_x_xdot_z;
external_function_casadi * impl_dae_jac_x_xdot_u_z;
{% else: %}
external_function_param_casadi * impl_dae_fun;
external_function_param_casadi * impl_dae_fun_jac_x_xdot_z;
external_function_param_casadi * impl_dae_jac_x_xdot_u_z;
{% endif %}
{% endif %}

int acados_create();
int acados_solve();
int acados_free();

ocp_nlp_in * acados_get_nlp_in();
ocp_nlp_out * acados_get_nlp_out();
ocp_nlp_solver * acados_get_nlp_solver();
ocp_nlp_solver_config * acados_get_nlp_config();
void * acados_get_nlp_opts();
ocp_nlp_dims * acados_get_nlp_dims();

#endif  // ACADOS_SOLVER_{{ra.model_name.upper()}}_H_
