
#ifndef EXAMPLES_C_CRANE_DAE
#define EXAMPLES_C_CRANE_DAE

#ifdef __cplusplus
extern "C" {
#endif

// this is a crane model with an artificially added algebraic equation to test gnsf & dae integrators

/* implicit ODE */

// implicit ODE
int crane_dae_impl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int crane_dae_impl_ode_fun_work(int *, int *, int *, int *);
const int *crane_dae_impl_ode_fun_sparsity_in(int);
const int *crane_dae_impl_ode_fun_sparsity_out(int);
int crane_dae_impl_ode_fun_n_in();
int crane_dae_impl_ode_fun_n_out();

// implicit ODE
int crane_dae_impl_ode_fun_jac_x_xdot(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int crane_dae_impl_ode_fun_jac_x_xdot_work(int *, int *, int *, int *);
const int *crane_dae_impl_ode_fun_jac_x_xdot_sparsity_in(int);
const int *crane_dae_impl_ode_fun_jac_x_xdot_sparsity_out(int);
int crane_dae_impl_ode_fun_jac_x_xdot_n_in();
int crane_dae_impl_ode_fun_jac_x_xdot_n_out();

// implicit ODE
int crane_dae_impl_ode_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int crane_dae_impl_ode_jac_x_xdot_u_work(int *, int *, int *, int *);
const int *crane_dae_impl_ode_jac_x_xdot_u_sparsity_in(int);
const int *crane_dae_impl_ode_jac_x_xdot_u_sparsity_out(int);
int crane_dae_impl_ode_jac_x_xdot_u_n_in();
int crane_dae_impl_ode_jac_x_xdot_u_n_out();

// implicit ODE - for lifted_irk
int crane_dae_impl_ode_fun_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int crane_dae_impl_ode_fun_jac_x_xdot_u_work(int *, int *, int *, int *);
const int *crane_dae_impl_ode_fun_jac_x_xdot_u_sparsity_in(int);
const int *crane_dae_impl_ode_fun_jac_x_xdot_u_sparsity_out(int);
int crane_dae_impl_ode_fun_jac_x_xdot_u_n_in();
int crane_dae_impl_ode_fun_jac_x_xdot_u_n_out();

/* GNSF Functions */
// used to import model matrices
int        crane_dae_get_matrices_fun(const double** arg, double** res, int* iw, double* w, void *mem);
int        crane_dae_get_matrices_fun_work(int *, int *, int *, int *);
const int *crane_dae_get_matrices_fun_sparsity_in(int);
const int *crane_dae_get_matrices_fun_sparsity_out(int);
int        crane_dae_get_matrices_fun_n_in();
int        crane_dae_get_matrices_fun_n_out();

// phi_fun
int        crane_dae_phi_fun(const double** arg, double** res, int* iw, double* w, void *mem);
int        crane_dae_phi_fun_work(int *, int *, int *, int *);
const int *crane_dae_phi_fun_sparsity_in(int);
const int *crane_dae_phi_fun_sparsity_out(int);
int        crane_dae_phi_fun_n_in();
int        crane_dae_phi_fun_n_out();

// phi_fun_jac_y
int        crane_dae_phi_fun_jac_y(const double** arg, double** res, int* iw, double* w, void *mem);
int        crane_dae_phi_fun_jac_y_work(int *, int *, int *, int *);
const int *crane_dae_phi_fun_jac_y_sparsity_in(int);
const int *crane_dae_phi_fun_jac_y_sparsity_out(int);
int        crane_dae_phi_fun_jac_y_n_in();
int        crane_dae_phi_fun_jac_y_n_out();

// phi_jac_y_uhat
int        crane_dae_phi_jac_y_uhat(const double** arg, double** res, int* iw, double* w, void *mem);
int        crane_dae_phi_jac_y_uhat_work(int *, int *, int *, int *);
const int *crane_dae_phi_jac_y_uhat_sparsity_in(int);
const int *crane_dae_phi_jac_y_uhat_sparsity_out(int);
int        crane_dae_phi_jac_y_uhat_n_in();
int        crane_dae_phi_jac_y_uhat_n_out();

// f_lo_fun_jac_x1k1uz
int        crane_dae_f_lo_fun_jac_x1k1uz(const double** arg, double** res, int* iw, double* w, void *mem);
int        crane_dae_f_lo_fun_jac_x1k1uz_work(int *, int *, int *, int *);
const int *crane_dae_f_lo_fun_jac_x1k1uz_sparsity_in(int);
const int *crane_dae_f_lo_fun_jac_x1k1uz_sparsity_out(int);
int        crane_dae_f_lo_fun_jac_x1k1uz_n_in();
int        crane_dae_f_lo_fun_jac_x1k1uz_n_out();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_CRANE_DAE
