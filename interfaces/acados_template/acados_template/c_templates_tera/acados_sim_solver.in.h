#ifndef ACADOS_SIM_{{model.name}}_H_
#define ACADOS_SIM_{{model.name}}_H_

#include "acados_c/sim_interface.h"
#include "acados_c/external_function_interface.h"

#ifdef __cplusplus
extern "C" {
#endif

int {{model.name}}_acados_sim_create();
int {{model.name}}_acados_sim_solve();
int {{model.name}}_acados_sim_free();

sim_config  * {{model.name}}_acados_get_sim_config();
sim_in      * {{model.name}}_acados_get_sim_in();
sim_out     * {{model.name}}_acados_get_sim_out();
void        * {{model.name}}_acados_get_sim_dims();
sim_opts    * {{model.name}}_acados_get_sim_opts();
sim_solver  * {{model.name}}_acados_get_sim_solver();

// ** global data **
extern sim_config  * {{model.name}}_sim_config;
extern sim_in      * {{model.name}}_sim_in;
extern sim_out     * {{model.name}}_sim_out; 
extern void        * {{model.name}}_sim_dims;
extern sim_opts    * {{model.name}}_sim_opts;
extern sim_solver  * {{model.name}}_sim_solver; 

#ifdef __cplusplus
}
#endif

{% if solver_config.integrator_type == "ERK" %}
extern external_function_param_casadi * sim_forw_vde_casadi;
extern external_function_param_casadi * sim_expl_ode_fun_casadi;
{% if solver_config.hessian_approx == "EXACT" %}
extern external_function_param_casadi * sim_hess_vde_casadi;
{% endif %}
{% else %}
{% if solver_config.integrator_type == "IRK" %}
extern external_function_param_casadi * sim_impl_dae_fun;
extern external_function_param_casadi * sim_impl_dae_fun_jac_x_xdot_z;
extern external_function_param_casadi * sim_impl_dae_jac_x_xdot_u_z;
{% endif %}
{% endif %}

#endif  // ACADOS_SIM_{{model.name}}_H_
