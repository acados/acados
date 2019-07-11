// TODO templetize with casadi function names !!!

// explicit ode function
int sim_model_dyn_expl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int sim_model_dyn_expl_ode_fun_work(int *, int *, int *, int *);
const int *sim_model_dyn_expl_ode_fun_sparsity_in(int);
const int *sim_model_dyn_expl_ode_fun_sparsity_out(int);
int sim_model_dyn_expl_ode_fun_n_in();
int sim_model_dyn_expl_ode_fun_n_out();

// explicit vde forward
int sim_model_dyn_expl_vde_for(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int sim_model_dyn_expl_vde_for_work(int *, int *, int *, int *);
const int *sim_model_dyn_expl_vde_for_sparsity_in(int);
const int *sim_model_dyn_expl_vde_for_sparsity_out(int);
int sim_model_dyn_expl_vde_for_n_in();
int sim_model_dyn_expl_vde_for_n_out();

// explicit vde adjoint
int sim_model_dyn_expl_vde_adj(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int sim_model_dyn_expl_vde_adj_work(int *, int *, int *, int *);
const int *sim_model_dyn_expl_vde_adj_sparsity_in(int);
const int *sim_model_dyn_expl_vde_adj_sparsity_out(int);
int sim_model_dyn_expl_vde_adj_n_in();
int sim_model_dyn_expl_vde_adj_n_out();

// explicit ode hessian
int sim_model_dyn_expl_ode_hes(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int sim_model_dyn_expl_ode_hes_work(int *, int *, int *, int *);
const int *sim_model_dyn_expl_ode_hes_sparsity_in(int);
const int *sim_model_dyn_expl_ode_hes_sparsity_out(int);
int sim_model_dyn_expl_ode_hes_n_in();
int sim_model_dyn_expl_ode_hes_n_out();

// implicit ode function
int sim_model_dyn_impl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int sim_model_dyn_impl_ode_fun_work(int *, int *, int *, int *);
const int *sim_model_dyn_impl_ode_fun_sparsity_in(int);
const int *sim_model_dyn_impl_ode_fun_sparsity_out(int);
int sim_model_dyn_impl_ode_fun_n_in();
int sim_model_dyn_impl_ode_fun_n_out();

// implicit ode function jacobian_x_xdot
int sim_model_dyn_impl_ode_fun_jac_x_xdot(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int sim_model_dyn_impl_ode_fun_jac_x_xdot_work(int *, int *, int *, int *);
const int *sim_model_dyn_impl_ode_fun_jac_x_xdot_sparsity_in(int);
const int *sim_model_dyn_impl_ode_fun_jac_x_xdot_sparsity_out(int);
int sim_model_dyn_impl_ode_fun_jac_x_xdot_n_in();
int sim_model_dyn_impl_ode_fun_jac_x_xdot_n_out();

// implicit ode jacobian_x_xdot_u
int sim_model_dyn_impl_ode_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int sim_model_dyn_impl_ode_jac_x_xdot_u_work(int *, int *, int *, int *);
const int *sim_model_dyn_impl_ode_jac_x_xdot_u_sparsity_in(int);
const int *sim_model_dyn_impl_ode_jac_x_xdot_u_sparsity_out(int);
int sim_model_dyn_impl_ode_jac_x_xdot_u_n_in();
int sim_model_dyn_impl_ode_jac_x_xdot_u_n_out();

// implicit ode hessian
int sim_model_dyn_impl_ode_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int sim_model_dyn_impl_ode_hess_work(int *, int *, int *, int *);
const int *sim_model_dyn_impl_ode_hess_sparsity_in(int);
const int *sim_model_dyn_impl_ode_hess_sparsity_out(int);
int sim_model_dyn_impl_ode_hess_n_in();
int sim_model_dyn_impl_ode_hess_n_out();

// gnsf
int sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz_work(int *, int *, int *, int *);
const int *sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz_sparsity_in(int);
const int *sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz_sparsity_out(int);
int sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz_n_in();
int sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz_n_out();

// gnsf
int sim_model_dyn_gnsf_get_matrices_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int sim_model_dyn_gnsf_get_matrices_fun_work(int *, int *, int *, int *);
const int *sim_model_dyn_gnsf_get_matrices_fun_sparsity_in(int);
const int *sim_model_dyn_gnsf_get_matrices_fun_sparsity_out(int);
int sim_model_dyn_gnsf_get_matrices_fun_n_in();
int sim_model_dyn_gnsf_get_matrices_fun_n_out();

// gnsf
int sim_model_dyn_gnsf_phi_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int sim_model_dyn_gnsf_phi_fun_work(int *, int *, int *, int *);
const int *sim_model_dyn_gnsf_phi_fun_sparsity_in(int);
const int *sim_model_dyn_gnsf_phi_fun_sparsity_out(int);
int sim_model_dyn_gnsf_phi_fun_n_in();
int sim_model_dyn_gnsf_phi_fun_n_out();

// gnsf
int sim_model_dyn_gnsf_phi_fun_jac_y(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int sim_model_dyn_gnsf_phi_fun_jac_y_work(int *, int *, int *, int *);
const int *sim_model_dyn_gnsf_phi_fun_jac_y_sparsity_in(int);
const int *sim_model_dyn_gnsf_phi_fun_jac_y_sparsity_out(int);
int sim_model_dyn_gnsf_phi_fun_jac_y_n_in();
int sim_model_dyn_gnsf_phi_fun_jac_y_n_out();

// gnsf
int sim_model_dyn_gnsf_phi_jac_y_uhat(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int sim_model_dyn_gnsf_phi_jac_y_uhat_work(int *, int *, int *, int *);
const int *sim_model_dyn_gnsf_phi_jac_y_uhat_sparsity_in(int);
const int *sim_model_dyn_gnsf_phi_jac_y_uhat_sparsity_out(int);
int sim_model_dyn_gnsf_phi_jac_y_uhat_n_in();
int sim_model_dyn_gnsf_phi_jac_y_uhat_n_out();

