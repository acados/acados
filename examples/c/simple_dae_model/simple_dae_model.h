
#ifndef EXAMPLES_C_SIMPLE_DAE
#define EXAMPLES_C_SIMPLE_DAE

#ifdef __cplusplus
extern "C" {
#endif

// implicit ODE
int simple_dae_impl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int simple_dae_impl_ode_fun_work(int *, int *, int *, int *);
const int *simple_dae_impl_ode_fun_sparsity_in(int);
const int *simple_dae_impl_ode_fun_sparsity_out(int);
int simple_dae_impl_ode_fun_n_in();
int simple_dae_impl_ode_fun_n_out();

// implicit ODE
int simple_dae_impl_ode_fun_jac_x_xdot_z(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int simple_dae_impl_ode_fun_jac_x_xdot_z_work(int *, int *, int *, int *);
const int *simple_dae_impl_ode_fun_jac_x_xdot_z_sparsity_in(int);
const int *simple_dae_impl_ode_fun_jac_x_xdot_z_sparsity_out(int);
int simple_dae_impl_ode_fun_jac_x_xdot_z_n_in();
int simple_dae_impl_ode_fun_jac_x_xdot_z_n_out();

// implicit ODE
int simple_dae_impl_ode_jac_x_xdot_u_z(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int simple_dae_impl_ode_jac_x_xdot_u_z_work(int *, int *, int *, int *);
const int *simple_dae_impl_ode_jac_x_xdot_u_z_sparsity_in(int);
const int *simple_dae_impl_ode_jac_x_xdot_u_z_sparsity_out(int);
int simple_dae_impl_ode_jac_x_xdot_u_z_n_in();
int simple_dae_impl_ode_jac_x_xdot_u_z_n_out();

// implicit ODE - for new_lifted_irk
int simple_dae_impl_ode_fun_jac_x_xdot_u_z(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int simple_dae_impl_ode_fun_jac_x_xdot_u_z_work(int *, int *, int *, int *);
const int *simple_dae_impl_ode_fun_jac_x_xdot_u_z_sparsity_in(int);
const int *simple_dae_impl_ode_fun_jac_x_xdot_u_z_sparsity_out(int);
int simple_dae_impl_ode_fun_jac_x_xdot_u_z_n_in();
int simple_dae_impl_ode_fun_jac_x_xdot_u_z_n_out();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_SIMPLE_DAE
