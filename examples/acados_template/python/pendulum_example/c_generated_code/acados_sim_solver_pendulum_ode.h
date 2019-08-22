#ifndef ACADOS_SIM_pendulum_ode_H_
#define ACADOS_SIM_pendulum_ode_H_

#include "acados_c/sim_interface.h"
#include "acados_c/external_function_interface.h"

#ifdef __cplusplus
extern "C" {
#endif

int pendulum_ode_acados_sim_create();
int pendulum_ode_acados_sim_solve();
int pendulum_ode_acados_sim_free();

sim_config  * pendulum_ode_acados_get_sim_config();
sim_in      * pendulum_ode_acados_get_sim_in();
sim_out     * pendulum_ode_acados_get_sim_out();
void        * pendulum_ode_acados_get_sim_dims();
sim_opts    * pendulum_ode_acados_get_sim_opts();
sim_solver  * pendulum_ode_acados_get_sim_solver();

// ** global data **
//
extern sim_config  * pendulum_ode_sim_config;
extern sim_in      * pendulum_ode_sim_in;
extern sim_out     * pendulum_ode_sim_out; 
extern void        * pendulum_ode_sim_dims;
extern sim_opts    * pendulum_ode_sim_opts;
extern sim_solver  * pendulum_ode_sim_solver; 

#ifdef __cplusplus
}
#endif



extern external_function_casadi * sim_forw_vde_casadi;
extern external_function_casadi * sim_expl_ode_fun_casadi;




#endif  // ACADOS_SIM_pendulum_ode_H_