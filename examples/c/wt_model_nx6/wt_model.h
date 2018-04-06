#ifndef EXAMPLES_C_WT_MODEL_NX6_H_
#define EXAMPLES_C_WT_MODEL_NX6_H_

#ifdef __cplusplus
extern "C" {
#endif


/* explicit ODE */

// explicit ODE
int casadi_expl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int casadi_expl_ode_fun_work(int *, int *, int *, int *);
const int *casadi_expl_ode_fun_sparsity_in(int);
const int *casadi_expl_ode_fun_sparsity_out(int);
int casadi_expl_ode_fun_n_in();
int casadi_expl_ode_fun_n_out();

// explicit forward VDE
int casadi_expl_vde_for(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int casadi_expl_vde_for_work(int *, int *, int *, int *);
const int *casadi_expl_vde_for_sparsity_in(int);
const int *casadi_expl_vde_for_sparsity_out(int);
int casadi_expl_vde_for_n_in();
int casadi_expl_vde_for_n_out();


/* implicit ODE */

// implicit ODE
int casadi_impl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int casadi_impl_ode_fun_work(int *, int *, int *, int *);
const int *casadi_impl_ode_fun_sparsity_in(int);
const int *casadi_impl_ode_fun_sparsity_out(int);
int casadi_impl_ode_fun_n_in();
int casadi_impl_ode_fun_n_out();

// implicit ODE
int casadi_impl_ode_jac_x(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int casadi_impl_ode_jac_x_work(int *, int *, int *, int *);
const int *casadi_impl_ode_jac_x_sparsity_in(int);
const int *casadi_impl_ode_jac_x_sparsity_out(int);
int casadi_impl_ode_jac_x_n_in();
int casadi_impl_ode_jac_x_n_out();

// implicit ODE
int casadi_impl_ode_jac_xdot(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int casadi_impl_ode_jac_xdot_work(int *, int *, int *, int *);
const int *casadi_impl_ode_jac_xdot_sparsity_in(int);
const int *casadi_impl_ode_jac_xdot_sparsity_out(int);
int casadi_impl_ode_jac_xdot_n_in();
int casadi_impl_ode_jac_xdot_n_out();

// implicit ODE
int casadi_impl_ode_jac_u(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int casadi_impl_ode_jac_u_work(int *, int *, int *, int *);
const int *casadi_impl_ode_jac_u_sparsity_in(int);
const int *casadi_impl_ode_jac_u_sparsity_out(int);
int casadi_impl_ode_jac_u_n_in();
int casadi_impl_ode_jac_u_n_out();

// implicit ODE
int casadi_impl_ode_fun_jac_x_xdot(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int casadi_impl_ode_fun_jac_x_xdot_work(int *, int *, int *, int *);
const int *casadi_impl_ode_fun_jac_x_xdot_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_sparsity_out(int);
int casadi_impl_ode_fun_jac_x_xdot_n_in();
int casadi_impl_ode_fun_jac_x_xdot_n_out();

// implicit ODE
int casadi_impl_ode_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int casadi_impl_ode_jac_x_xdot_u_work(int *, int *, int *, int *);
const int *casadi_impl_ode_jac_x_xdot_u_sparsity_in(int);
const int *casadi_impl_ode_jac_x_xdot_u_sparsity_out(int);
int casadi_impl_ode_jac_x_xdot_u_n_in();
int casadi_impl_ode_jac_x_xdot_u_n_out();

// implicit ODE
int casadi_impl_ode_jac_x_u(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int casadi_impl_ode_jac_x_u_work(int *, int *, int *, int *);
const int *casadi_impl_ode_jac_x_u_sparsity_in(int);
const int *casadi_impl_ode_jac_x_u_sparsity_out(int);
int casadi_impl_ode_jac_x_u_n_in();
int casadi_impl_ode_jac_x_u_n_out();

/* GNSF Functions */
// used to import integers & double matrices
int get_ints_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int get_matrices_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
// int But_KK_YY_ZZ_LO_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

// Phi_inc_dy_fun
int casadi_phi_fun_jac_y(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int casadi_phi_fun_jac_y_work(int *, int *, int *, int *);
const int *casadi_phi_fun_jac_y_sparsity_in(int);
const int *casadi_phi_fun_jac_y_sparsity_out(int);
int casadi_phi_fun_jac_y_n_in();
int casadi_phi_fun_jac_y_n_out();


//Phi_jac_y_fun
int casadi_phi_jac_y(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int casadi_phi_jac_y_work(int *, int *, int *, int *);
const int *casadi_phi_jac_y_sparsity_in(int);
const int *casadi_phi_jac_y_sparsity_out(int);
int casadi_phi_jac_y_n_in();
int casadi_phi_jac_y_n_out();

// f_LO_inc_J_x1k1uz_fun
int casadi_f_LO_inc_J_x1k1uz_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int casadi_f_LO_inc_J_x1k1uz_fun_work(int *, int *, int *, int *);
const int *casadi_f_LO_inc_J_x1k1uz_fun_sparsity_in(int);
const int *casadi_f_LO_inc_J_x1k1uz_fun_sparsity_out(int);
int casadi_f_LO_inc_J_x1k1uz_fun_n_in();
int casadi_f_LO_inc_J_x1k1uz_fun_n_out();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_WT_MODEL_NX6_H_
