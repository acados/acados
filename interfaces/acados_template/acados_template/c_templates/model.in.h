#ifndef {{ ra.model_name.upper() }}_MODEL
#define {{ ra.model_name.upper() }}_MODEL

#ifdef __cplusplus
extern "C" {
#endif

{% if ra.solver_config.integrator_type != 'ERK' %}
// implicit ODE
int {{ ra.model_name }}_impl_dae_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ra.model_name }}_impl_dae_fun_work(int *, int *, int *, int *);
const int *{{ ra.model_name }}_impl_dae_fun_sparsity_in(int);
const int *{{ ra.model_name }}_impl_dae_fun_sparsity_out(int);
int {{ ra.model_name }}_impl_dae_fun_n_in();
int {{ ra.model_name }}_impl_dae_fun_n_out();

// implicit ODE
int {{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z_work(int *, int *, int *, int *);
const int *{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z_sparsity_in(int);
const int *{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z_sparsity_out(int);
int {{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z_n_in();
int {{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z_n_out();

// implicit ODE
int {{ ra.model_name }}_impl_dae_jac_x_xdot_u_z(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ra.model_name }}_impl_dae_jac_x_xdot_u_z_work(int *, int *, int *, int *);
const int *{{ ra.model_name }}_impl_dae_jac_x_xdot_u_z_sparsity_in(int);
const int *{{ ra.model_name }}_impl_dae_jac_x_xdot_u_z_sparsity_out(int);
int {{ ra.model_name }}_impl_dae_jac_x_xdot_u_z_n_in();
int {{ ra.model_name }}_impl_dae_jac_x_xdot_u_z_n_out();

// // implicit ODE - for lifted_irk
// int {{ ra.model_name }}_impl_dae_fun_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
// int {{ ra.model_name }}_impl_dae_fun_jac_x_xdot_u_work(int *, int *, int *, int *);
// const int *{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_u_sparsity_in(int);
// const int *{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_u_sparsity_out(int);
// int {{ ra.model_name }}_impl_dae_fun_jac_x_xdot_u_n_in();
// int {{ ra.model_name }}_impl_dae_fun_jac_x_xdot_u_n_out();

// // implicit ODE
// int {{ ra.model_name }}_impl_dae_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
// int {{ ra.model_name }}_impl_dae_hess_work(int *, int *, int *, int *);
// const int *{{ ra.model_name }}_impl_dae_hess_sparsity_in(int);
// const int *{{ ra.model_name }}_impl_dae_hess_sparsity_out(int);
// int {{ ra.model_name }}_impl_dae_hess_n_in();
// int {{ ra.model_name }}_impl_dae_hess_n_out();

// /* GNSF Functions */
// // used to import model matrices
// int        {{ ra.model_name }}_get_matrices_fun(const double** arg, double** res, int* iw, double* w, void *mem);
// int        {{ ra.model_name }}_get_matrices_fun_work(int *, int *, int *, int *);
// const int *{{ ra.model_name }}_get_matrices_fun_sparsity_in(int);
// const int *{{ ra.model_name }}_get_matrices_fun_sparsity_out(int);
// int        {{ ra.model_name }}_get_matrices_fun_n_in();
// int        {{ ra.model_name }}_get_matrices_fun_n_out();

// // phi_fun
// int        {{ ra.model_name }}_phi_fun(const double** arg, double** res, int* iw, double* w, void *mem);
// int        {{ ra.model_name }}_phi_fun_work(int *, int *, int *, int *);
// const int *{{ ra.model_name }}_phi_fun_sparsity_in(int);
// const int *{{ ra.model_name }}_phi_fun_sparsity_out(int);
// int        {{ ra.model_name }}_phi_fun_n_in();
// int        {{ ra.model_name }}_phi_fun_n_out();

// // phi_fun_jac_y
// int        {{ ra.model_name }}_phi_fun_jac_y(const double** arg, double** res, int* iw, double* w, void *mem);
// int        {{ ra.model_name }}_phi_fun_jac_y_work(int *, int *, int *, int *);
// const int *{{ ra.model_name }}_phi_fun_jac_y_sparsity_in(int);
// const int *{{ ra.model_name }}_phi_fun_jac_y_sparsity_out(int);
// int        {{ ra.model_name }}_phi_fun_jac_y_n_in();
// int        {{ ra.model_name }}_phi_fun_jac_y_n_out();

// // phi_jac_y_uhat
// int        {{ ra.model_name }}_phi_jac_y_uhat(const double** arg, double** res, int* iw, double* w, void *mem);
// int        {{ ra.model_name }}_phi_jac_y_uhat_work(int *, int *, int *, int *);
// const int *{{ ra.model_name }}_phi_jac_y_uhat_sparsity_in(int);
// const int *{{ ra.model_name }}_phi_jac_y_uhat_sparsity_out(int);
// int        {{ ra.model_name }}_phi_jac_y_uhat_n_in();
// int        {{ ra.model_name }}_phi_jac_y_uhat_n_out();

// // f_lo_fun_jac_x1k1uz
// int        {{ ra.model_name }}_f_lo_fun_jac_x1k1uz(const double** arg, double** res, int* iw, double* w, void *mem);
// int        {{ ra.model_name }}_f_lo_fun_jac_x1k1uz_work(int *, int *, int *, int *);
// const int *{{ ra.model_name }}_f_lo_fun_jac_x1k1uz_sparsity_in(int);
// const int *{{ ra.model_name }}_f_lo_fun_jac_x1k1uz_sparsity_out(int);
// int        {{ ra.model_name }}_f_lo_fun_jac_x1k1uz_n_in();
// int        {{ ra.model_name }}_f_lo_fun_jac_x1k1uz_n_out();

{% else %}
/* explicit ODE */

// explicit ODE
int {{ ra.model_name }}_expl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ra.model_name }}_expl_ode_fun_work(int *, int *, int *, int *);
const int *{{ ra.model_name }}_expl_ode_fun_sparsity_in(int);
const int *{{ ra.model_name }}_expl_ode_fun_sparsity_out(int);
int {{ ra.model_name }}_expl_ode_fun_n_in();
int {{ ra.model_name }}_expl_ode_fun_n_out();

// explicit forward VDE
int {{ ra.model_name }}_expl_vde_forw(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ra.model_name }}_expl_vde_forw_work(int *, int *, int *, int *);
const int *{{ ra.model_name }}_expl_vde_forw_sparsity_in(int);
const int *{{ ra.model_name }}_expl_vde_forw_sparsity_out(int);
int {{ ra.model_name }}_expl_vde_forw_n_in();
int {{ ra.model_name }}_expl_vde_forw_n_out();

// explicit adjoint VDE
int {{ ra.model_name }}_expl_vde_adj(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ra.model_name }}_expl_vde_adj_work(int *, int *, int *, int *);
const int *{{ ra.model_name }}_expl_vde_adj_sparsity_in(int);
const int *{{ ra.model_name }}_expl_vde_adj_sparsity_out(int);
int {{ ra.model_name }}_expl_vde_adj_n_in();
int {{ ra.model_name }}_expl_vde_adj_n_out();

// explicit adjoint ODE jac
int {{ ra.model_name }}_expl_ode_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ra.model_name }}_expl_ode_hess_work(int *, int *, int *, int *);
const int *{{ ra.model_name }}_expl_ode_hess_sparsity_in(int);
const int *{{ ra.model_name }}_expl_ode_hess_sparsity_out(int);
int {{ ra.model_name }}_expl_ode_hess_n_in();
int {{ ra.model_name }}_expl_ode_hess_n_out();

{% endif %}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // {{ ra.model_name.upper() }}_MODEL
