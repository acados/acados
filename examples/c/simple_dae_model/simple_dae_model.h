
#ifndef EXAMPLES_C_SIMPLE_DAE
#define EXAMPLES_C_SIMPLE_DAE

#ifdef __cplusplus
extern "C" {
#endif

// implicit ODE
int casadi_impl_ode_fun_simple_dae(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_simple_dae_work(int *, int *, int *, int *);
const int *casadi_impl_ode_fun_simple_dae_sparsity_in(int);
const int *casadi_impl_ode_fun_simple_dae_sparsity_out(int);
int casadi_impl_ode_fun_simple_dae_n_in();
int casadi_impl_ode_fun_simple_dae_n_out();

// implicit ODE
int casadi_impl_ode_fun_jac_x_xdot_z_simple_dae(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_z_simple_dae_work(int *, int *, int *, int *);
const int *casadi_impl_ode_fun_jac_x_xdot_z_simple_dae_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_z_simple_dae_sparsity_out(int);
int casadi_impl_ode_fun_jac_x_xdot_z_simple_dae_n_in();
int casadi_impl_ode_fun_jac_x_xdot_z_simple_dae_n_out();

// implicit ODE
int casadi_impl_ode_jac_x_xdot_u_z_simple_dae(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_jac_x_xdot_u_z_simple_dae_work(int *, int *, int *, int *);
const int *casadi_impl_ode_jac_x_xdot_u_z_simple_dae_sparsity_in(int);
const int *casadi_impl_ode_jac_x_xdot_u_z_simple_dae_sparsity_out(int);
int casadi_impl_ode_jac_x_xdot_u_z_simple_dae_n_in();
int casadi_impl_ode_jac_x_xdot_u_z_simple_dae_n_out();

// implicit ODE - for new_lifted_irk
int casadi_impl_ode_fun_jac_x_xdot_u_z_simple_dae(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_u_z_simple_dae_work(int *, int *, int *, int *);
const int *casadi_impl_ode_fun_jac_x_xdot_u_z_simple_dae_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_u_z_simple_dae_sparsity_out(int);
int casadi_impl_ode_fun_jac_x_xdot_u_z_simple_dae_n_in();
int casadi_impl_ode_fun_jac_x_xdot_u_z_simple_dae_n_out();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_SIMPLE_DAE
