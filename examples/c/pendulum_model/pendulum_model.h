
#ifndef EXAMPLES_C_PENDULUM_MODEL_H_
#define EXAMPLES_C_PENDULUM_MODEL_H_

#ifdef __cplusplus
extern "C" {
#endif


/* explicit ODE */

// explicit ODE
int pendulum_ode_expl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_expl_ode_fun_work(int *, int *, int *, int *);
const int *pendulum_ode_expl_ode_fun_sparsity_in(int);
const int *pendulum_ode_expl_ode_fun_sparsity_out(int);
int pendulum_ode_expl_ode_fun_n_in();
int pendulum_ode_expl_ode_fun_n_out();

// explicit forward VDE
int pendulum_ode_expl_vde_forw(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_expl_vde_forw_work(int *, int *, int *, int *);
const int *pendulum_ode_expl_vde_forw_sparsity_in(int);
const int *pendulum_ode_expl_vde_forw_sparsity_out(int);
int pendulum_ode_expl_vde_forw_n_in();
int pendulum_ode_expl_vde_forw_n_out();

// explicit adjoint VDE
int pendulum_ode_expl_vde_adj(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_expl_vde_adj_work(int *, int *, int *, int *);
const int *pendulum_ode_expl_vde_adj_sparsity_in(int);
const int *pendulum_ode_expl_vde_adj_sparsity_out(int);
int pendulum_ode_expl_vde_adj_n_in();
int pendulum_ode_expl_vde_adj_n_out();

// explicit adjoint ODE jac
int pendulum_ode_expl_ode_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_expl_ode_hess_work(int *, int *, int *, int *);
const int *pendulum_ode_expl_ode_hess_sparsity_in(int);
const int *pendulum_ode_expl_ode_hess_sparsity_out(int);
int pendulum_ode_expl_ode_hess_n_in();
int pendulum_ode_expl_ode_hess_n_out();


/* implicit ODE */

// implicit ODE
int pendulum_ode_impl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_impl_ode_fun_work(int *, int *, int *, int *);
const int *pendulum_ode_impl_ode_fun_sparsity_in(int);
const int *pendulum_ode_impl_ode_fun_sparsity_out(int);
int pendulum_ode_impl_ode_fun_n_in();
int pendulum_ode_impl_ode_fun_n_out();

// implicit ODE
int pendulum_ode_impl_ode_fun_jac_x_xdot(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_impl_ode_fun_jac_x_xdot_work(int *, int *, int *, int *);
const int *pendulum_ode_impl_ode_fun_jac_x_xdot_sparsity_in(int);
const int *pendulum_ode_impl_ode_fun_jac_x_xdot_sparsity_out(int);
int pendulum_ode_impl_ode_fun_jac_x_xdot_n_in();
int pendulum_ode_impl_ode_fun_jac_x_xdot_n_out();

// implicit ODE
int pendulum_ode_impl_ode_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_impl_ode_jac_x_xdot_u_work(int *, int *, int *, int *);
const int *pendulum_ode_impl_ode_jac_x_xdot_u_sparsity_in(int);
const int *pendulum_ode_impl_ode_jac_x_xdot_u_sparsity_out(int);
int pendulum_ode_impl_ode_jac_x_xdot_u_n_in();
int pendulum_ode_impl_ode_jac_x_xdot_u_n_out();

// implicit ODE - for lifted_irk
int pendulum_ode_impl_ode_fun_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_impl_ode_fun_jac_x_xdot_u_work(int *, int *, int *, int *);
const int *pendulum_ode_impl_ode_fun_jac_x_xdot_u_sparsity_in(int);
const int *pendulum_ode_impl_ode_fun_jac_x_xdot_u_sparsity_out(int);
int pendulum_ode_impl_ode_fun_jac_x_xdot_u_n_in();
int pendulum_ode_impl_ode_fun_jac_x_xdot_u_n_out();

int pendulum_ode_impl_ode_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_impl_ode_hess_work(int *, int *, int *, int *);
const int *pendulum_ode_impl_ode_hess_sparsity_in(int);
const int *pendulum_ode_impl_ode_hess_sparsity_out(int);
int pendulum_ode_impl_ode_hess_n_in();
int pendulum_ode_impl_ode_hess_n_out();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_PENDULUM_MODEL_H_