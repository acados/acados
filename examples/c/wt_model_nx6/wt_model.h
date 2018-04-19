#ifndef EXAMPLES_C_WT_MODEL_NX6_H_
#define EXAMPLES_C_WT_MODEL_NX6_H_

#ifdef __cplusplus
extern "C" {
#endif


/* explicit ODE */

// explicit ODE
int casadi_expl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_expl_ode_fun_work(int *, int *, int *, int *);
const int *casadi_expl_ode_fun_sparsity_in(int);
const int *casadi_expl_ode_fun_sparsity_out(int);
int casadi_expl_ode_fun_n_in();
int casadi_expl_ode_fun_n_out();

// explicit forward VDE
int casadi_expl_vde_for(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_expl_vde_for_work(int *, int *, int *, int *);
const int *casadi_expl_vde_for_sparsity_in(int);
const int *casadi_expl_vde_for_sparsity_out(int);
int casadi_expl_vde_for_n_in();
int casadi_expl_vde_for_n_out();

// explicit adjoint VDE
int casadi_expl_vde_adj(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_expl_vde_adj_work(int *, int *, int *, int *);
const int *casadi_expl_vde_adj_sparsity_in(int);
const int *casadi_expl_vde_adj_sparsity_out(int);
int casadi_expl_vde_adj_n_in();
int casadi_expl_vde_adj_n_out();


/* implicit ODE */

// implicit ODE
int casadi_impl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_work(int *, int *, int *, int *);
const int *casadi_impl_ode_fun_sparsity_in(int);
const int *casadi_impl_ode_fun_sparsity_out(int);
int casadi_impl_ode_fun_n_in();
int casadi_impl_ode_fun_n_out();

// implicit ODE
int casadi_impl_ode_jac_x(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_jac_x_work(int *, int *, int *, int *);
const int *casadi_impl_ode_jac_x_sparsity_in(int);
const int *casadi_impl_ode_jac_x_sparsity_out(int);
int casadi_impl_ode_jac_x_n_in();
int casadi_impl_ode_jac_x_n_out();

// implicit ODE
int casadi_impl_ode_jac_xdot(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_jac_xdot_work(int *, int *, int *, int *);
const int *casadi_impl_ode_jac_xdot_sparsity_in(int);
const int *casadi_impl_ode_jac_xdot_sparsity_out(int);
int casadi_impl_ode_jac_xdot_n_in();
int casadi_impl_ode_jac_xdot_n_out();

// implicit ODE
int casadi_impl_ode_jac_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_jac_u_work(int *, int *, int *, int *);
const int *casadi_impl_ode_jac_u_sparsity_in(int);
const int *casadi_impl_ode_jac_u_sparsity_out(int);
int casadi_impl_ode_jac_u_n_in();
int casadi_impl_ode_jac_u_n_out();

// implicit ODE
int casadi_impl_ode_fun_jac_x_xdot(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_work(int *, int *, int *, int *);
const int *casadi_impl_ode_fun_jac_x_xdot_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_sparsity_out(int);
int casadi_impl_ode_fun_jac_x_xdot_n_in();
int casadi_impl_ode_fun_jac_x_xdot_n_out();

// implicit ODE
int casadi_impl_ode_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_jac_x_xdot_u_work(int *, int *, int *, int *);
const int *casadi_impl_ode_jac_x_xdot_u_sparsity_in(int);
const int *casadi_impl_ode_jac_x_xdot_u_sparsity_out(int);
int casadi_impl_ode_jac_x_xdot_u_n_in();
int casadi_impl_ode_jac_x_xdot_u_n_out();

// implicit ODE
int casadi_impl_ode_jac_x_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_jac_x_u_work(int *, int *, int *, int *);
const int *casadi_impl_ode_jac_x_u_sparsity_in(int);
const int *casadi_impl_ode_jac_x_u_sparsity_out(int);
int casadi_impl_ode_jac_x_u_n_in();
int casadi_impl_ode_jac_x_u_n_out();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_WT_MODEL_NX6_H_
