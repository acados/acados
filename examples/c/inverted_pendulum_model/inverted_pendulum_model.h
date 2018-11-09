
#ifndef EXAMPLES_C_INV_PENDULUM
#define EXAMPLES_C_INV_PENDULUM

#ifdef __cplusplus
extern "C" {
#endif

// this is a crane model with an artificially added algebraic equation to test gnsf & dae integrators
/* implicit ODE */

// implicit ODE
int inv_pendulum_impl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int inv_pendulum_impl_ode_fun_work(int *, int *, int *, int *);
const int *inv_pendulum_impl_ode_fun_sparsity_in(int);
const int *inv_pendulum_impl_ode_fun_sparsity_out(int);
int inv_pendulum_impl_ode_fun_n_in();
int inv_pendulum_impl_ode_fun_n_out();

// implicit ODE
int inv_pendulum_impl_ode_fun_jac_x_xdot(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int inv_pendulum_impl_ode_fun_jac_x_xdot_work(int *, int *, int *, int *);
const int *inv_pendulum_impl_ode_fun_jac_x_xdot_sparsity_in(int);
const int *inv_pendulum_impl_ode_fun_jac_x_xdot_sparsity_out(int);
int inv_pendulum_impl_ode_fun_jac_x_xdot_n_in();
int inv_pendulum_impl_ode_fun_jac_x_xdot_n_out();

// implicit ODE
int inv_pendulum_impl_ode_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int inv_pendulum_impl_ode_jac_x_xdot_u_work(int *, int *, int *, int *);
const int *inv_pendulum_impl_ode_jac_x_xdot_u_sparsity_in(int);
const int *inv_pendulum_impl_ode_jac_x_xdot_u_sparsity_out(int);
int inv_pendulum_impl_ode_jac_x_xdot_u_n_in();
int inv_pendulum_impl_ode_jac_x_xdot_u_n_out();

// implicit ODE - for lifted_irk
int inv_pendulum_impl_ode_fun_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int inv_pendulum_impl_ode_fun_jac_x_xdot_u_work(int *, int *, int *, int *);
const int *inv_pendulum_impl_ode_fun_jac_x_xdot_u_sparsity_in(int);
const int *inv_pendulum_impl_ode_fun_jac_x_xdot_u_sparsity_out(int);
int inv_pendulum_impl_ode_fun_jac_x_xdot_u_n_in();
int inv_pendulum_impl_ode_fun_jac_x_xdot_u_n_out();

// implicit ODE
int inv_pendulum_impl_ode_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int inv_pendulum_impl_ode_hess_work(int *, int *, int *, int *);
const int *inv_pendulum_impl_ode_hess_sparsity_in(int);
const int *inv_pendulum_impl_ode_hess_sparsity_out(int);
int inv_pendulum_impl_ode_hess_n_in();
int inv_pendulum_impl_ode_hess_n_out();

/* GNSF Functions */
// used to import model matrices
int        inv_pendulum_get_matrices_fun(const double** arg, double** res, int* iw, double* w, void *mem);
int        inv_pendulum_get_matrices_fun_work(int *, int *, int *, int *);
const int *inv_pendulum_get_matrices_fun_sparsity_in(int);
const int *inv_pendulum_get_matrices_fun_sparsity_out(int);
int        inv_pendulum_get_matrices_fun_n_in();
int        inv_pendulum_get_matrices_fun_n_out();

// phi_fun
int        inv_pendulum_phi_fun(const double** arg, double** res, int* iw, double* w, void *mem);
int        inv_pendulum_phi_fun_work(int *, int *, int *, int *);
const int *inv_pendulum_phi_fun_sparsity_in(int);
const int *inv_pendulum_phi_fun_sparsity_out(int);
int        inv_pendulum_phi_fun_n_in();
int        inv_pendulum_phi_fun_n_out();

// phi_fun_jac_y
int        inv_pendulum_phi_fun_jac_y(const double** arg, double** res, int* iw, double* w, void *mem);
int        inv_pendulum_phi_fun_jac_y_work(int *, int *, int *, int *);
const int *inv_pendulum_phi_fun_jac_y_sparsity_in(int);
const int *inv_pendulum_phi_fun_jac_y_sparsity_out(int);
int        inv_pendulum_phi_fun_jac_y_n_in();
int        inv_pendulum_phi_fun_jac_y_n_out();

// phi_jac_y_uhat
int        inv_pendulum_phi_jac_y_uhat(const double** arg, double** res, int* iw, double* w, void *mem);
int        inv_pendulum_phi_jac_y_uhat_work(int *, int *, int *, int *);
const int *inv_pendulum_phi_jac_y_uhat_sparsity_in(int);
const int *inv_pendulum_phi_jac_y_uhat_sparsity_out(int);
int        inv_pendulum_phi_jac_y_uhat_n_in();
int        inv_pendulum_phi_jac_y_uhat_n_out();

// f_lo_fun_jac_x1k1uz
int        inv_pendulum_f_lo_fun_jac_x1k1uz(const double** arg, double** res, int* iw, double* w, void *mem);
int        inv_pendulum_f_lo_fun_jac_x1k1uz_work(int *, int *, int *, int *);
const int *inv_pendulum_f_lo_fun_jac_x1k1uz_sparsity_in(int);
const int *inv_pendulum_f_lo_fun_jac_x1k1uz_sparsity_out(int);
int        inv_pendulum_f_lo_fun_jac_x1k1uz_n_in();
int        inv_pendulum_f_lo_fun_jac_x1k1uz_n_out();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_INV_PENDULUM
