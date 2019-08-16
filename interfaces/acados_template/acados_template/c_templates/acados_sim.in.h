#ifndef ACADOS_SOLVER_{{ocp.model_name}}_H_
#define ACADOS_SOLVER_{{ocp.model_name}}_H_

#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

#ifdef __cplusplus
extern "C" {
#endif

int acados_create();
int acados_solve();
int acados_free();

ocp_nlp_in * acados_get_nlp_in();
ocp_nlp_out * acados_get_nlp_out();
ocp_nlp_solver * acados_get_nlp_solver();
ocp_nlp_config * acados_get_nlp_config();
void * acados_get_nlp_opts();
ocp_nlp_dims * acados_get_nlp_dims();

#ifdef __cplusplus
}
#endif


// ** global data **
extern ocp_nlp_in * nlp_in;
extern ocp_nlp_out * nlp_out;
extern ocp_nlp_solver * nlp_solver;
extern void * nlp_opts;
extern ocp_nlp_plan * nlp_solver_plan;
extern ocp_nlp_config * nlp_config;
extern ocp_nlp_dims * nlp_dims;
{% if ocp.solver_config.integrator_type == "ERK" %}
{% if ocp.dims.np < 1 %}
extern external_function_casadi * forw_vde_casadi;
{% else %}
extern external_function_param_casadi * forw_vde_casadi;
{% endif %}
{% if ocp.solver_config.hessian_approx == "EXACT" %}
{% if ocp.dims.np < 1 %}
extern external_function_casadi * hess_vde_casadi;
{% else %}
extern external_function_param_casadi * hess_vde_casadi;
{% endif %}
{% endif %}
{% else %}
{% if ocp.solver_config.integrator_type == "IRK" %}
{% if ocp.dims.np < 1 %}
extern external_function_casadi * impl_dae_fun;
extern external_function_casadi * impl_dae_fun_jac_x_xdot_z;
extern external_function_casadi * impl_dae_jac_x_xdot_u_z;
{% else %}
extern external_function_param_casadi * impl_dae_fun;
extern external_function_param_casadi * impl_dae_fun_jac_x_xdot_z;
extern external_function_param_casadi * impl_dae_jac_x_xdot_u_z;
{% endif %}
{% endif %}
{% endif %}
{% if ocp.dims.npd > 0 %}
extern external_function_casadi * p_constraint;
{% endif %}
{% if ocp.dims.npd_e > 0 %}
extern external_function_casadi * p_constraint_e;
{% endif %}
{% if ocp.dims.nh > 0 %}
extern external_function_casadi * h_constraint;
{% endif %}
{% if ocp.dims.nh_e > 0 %}
extern external_function_casadi * h_constraint_e;
{% endif %}

#endif  // ACADOS_SOLVER_{{ocp.model_name}}_H_
