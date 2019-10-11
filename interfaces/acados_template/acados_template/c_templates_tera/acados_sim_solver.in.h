#ifndef ACADOS_SIM_{{ocp.model.name}}_H_
#define ACADOS_SIM_{{ocp.model.name}}_H_

#include "acados_c/sim_interface.h"
#include "acados_c/external_function_interface.h"

#ifdef __cplusplus
extern "C" {
#endif

int {{ocp.model.name}}_acados_sim_create();
int {{ocp.model.name}}_acados_sim_solve();
int {{ocp.model.name}}_acados_sim_free();

sim_config  * {{ocp.model.name}}_acados_get_sim_config();
sim_in      * {{ocp.model.name}}_acados_get_sim_in();
sim_out     * {{ocp.model.name}}_acados_get_sim_out();
void        * {{ocp.model.name}}_acados_get_sim_dims();
sim_opts    * {{ocp.model.name}}_acados_get_sim_opts();
sim_solver  * {{ocp.model.name}}_acados_get_sim_solver();

// ** global data **
extern sim_config  * {{ocp.model.name}}_sim_config;
extern sim_in      * {{ocp.model.name}}_sim_in;
extern sim_out     * {{ocp.model.name}}_sim_out; 
extern void        * {{ocp.model.name}}_sim_dims;
extern sim_opts    * {{ocp.model.name}}_sim_opts;
extern sim_solver  * {{ocp.model.name}}_sim_solver; 

#ifdef __cplusplus
}
#endif

{% if ocp.solver_config.integrator_type == "ERK" %}
{% if ocp.dims.np < 1 %}
extern external_function_casadi * sim_forw_vde_casadi;
extern external_function_casadi * sim_expl_ode_fun_casadi;
{% else %}
extern external_function_param_casadi * sim_forw_vde_casadi;
extern external_function_param_casadi * sim_expl_ode_fun_casadi;
{% endif %}
{% if ocp.solver_config.hessian_approx == "EXACT" %}
{% if ocp.dims.np < 1 %}
extern external_function_casadi * sim_hess_vde_casadi;
{% else %}
extern external_function_param_casadi * sim_hess_vde_casadi;
{% endif %}
{% endif %}
{% else %}
{% if ocp.solver_config.integrator_type == "IRK" %}
{% if ocp.dims.np < 1 %}
extern external_function_casadi * sim_impl_dae_fun;
extern external_function_casadi * sim_impl_dae_fun_jac_x_xdot_z;
extern external_function_casadi * sim_impl_dae_jac_x_xdot_u_z;
{% else %}
extern external_function_param_casadi * sim_impl_dae_fun;
extern external_function_param_casadi * sim_impl_dae_fun_jac_x_xdot_z;
extern external_function_param_casadi * sim_impl_dae_jac_x_xdot_u_z;
{% endif %}
{% endif %}
{% endif %}

#endif  // ACADOS_SIM_{{ocp.model.name}}_H_
