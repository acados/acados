/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */

#ifndef pendulum_ode_MODEL
#define pendulum_ode_MODEL

#ifdef __cplusplus
extern "C" {
#endif


  
// implicit ODE: function
int pendulum_ode_impl_dae_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_impl_dae_fun_work(int *, int *, int *, int *);
const int *pendulum_ode_impl_dae_fun_sparsity_in(int);
const int *pendulum_ode_impl_dae_fun_sparsity_out(int);
int pendulum_ode_impl_dae_fun_n_in(void);
int pendulum_ode_impl_dae_fun_n_out(void);

// implicit ODE: function + jacobians
int pendulum_ode_impl_dae_fun_jac_x_xdot_z(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_impl_dae_fun_jac_x_xdot_z_work(int *, int *, int *, int *);
const int *pendulum_ode_impl_dae_fun_jac_x_xdot_z_sparsity_in(int);
const int *pendulum_ode_impl_dae_fun_jac_x_xdot_z_sparsity_out(int);
int pendulum_ode_impl_dae_fun_jac_x_xdot_z_n_in(void);
int pendulum_ode_impl_dae_fun_jac_x_xdot_z_n_out(void);

// implicit ODE: jacobians only
int pendulum_ode_impl_dae_jac_x_xdot_u_z(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_impl_dae_jac_x_xdot_u_z_work(int *, int *, int *, int *);
const int *pendulum_ode_impl_dae_jac_x_xdot_u_z_sparsity_in(int);
const int *pendulum_ode_impl_dae_jac_x_xdot_u_z_sparsity_out(int);
int pendulum_ode_impl_dae_jac_x_xdot_u_z_n_in(void);
int pendulum_ode_impl_dae_jac_x_xdot_u_z_n_out(void);

// implicit ODE - for lifted_irk
int pendulum_ode_impl_dae_fun_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_impl_dae_fun_jac_x_xdot_u_work(int *, int *, int *, int *);
const int *pendulum_ode_impl_dae_fun_jac_x_xdot_u_sparsity_in(int);
const int *pendulum_ode_impl_dae_fun_jac_x_xdot_u_sparsity_out(int);
int pendulum_ode_impl_dae_fun_jac_x_xdot_u_n_in(void);
int pendulum_ode_impl_dae_fun_jac_x_xdot_u_n_out(void);
  



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // pendulum_ode_MODEL
