// TODO templetize with casadi function names !!!

// explicit ode function
int ocp_model_dyn_expl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_dyn_expl_ode_fun_work(int *, int *, int *, int *);
const int *ocp_model_dyn_expl_ode_fun_sparsity_in(int);
const int *ocp_model_dyn_expl_ode_fun_sparsity_out(int);
int ocp_model_dyn_expl_ode_fun_n_in();
int ocp_model_dyn_expl_ode_fun_n_out();

// explicit vde forward
int ocp_model_dyn_expl_vde_for(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_dyn_expl_vde_for_work(int *, int *, int *, int *);
const int *ocp_model_dyn_expl_vde_for_sparsity_in(int);
const int *ocp_model_dyn_expl_vde_for_sparsity_out(int);
int ocp_model_dyn_expl_vde_for_n_in();
int ocp_model_dyn_expl_vde_for_n_out();

// explicit vde adjoint
int ocp_model_dyn_expl_vde_adj(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_dyn_expl_vde_adj_work(int *, int *, int *, int *);
const int *ocp_model_dyn_expl_vde_adj_sparsity_in(int);
const int *ocp_model_dyn_expl_vde_adj_sparsity_out(int);
int ocp_model_dyn_expl_vde_adj_n_in();
int ocp_model_dyn_expl_vde_adj_n_out();

// explicit ode hessian
int ocp_model_dyn_expl_ode_hes(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_dyn_expl_ode_hes_work(int *, int *, int *, int *);
const int *ocp_model_dyn_expl_ode_hes_sparsity_in(int);
const int *ocp_model_dyn_expl_ode_hes_sparsity_out(int);
int ocp_model_dyn_expl_ode_hes_n_in();
int ocp_model_dyn_expl_ode_hes_n_out();

// implicit ode function
int ocp_model_dyn_impl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_dyn_impl_ode_fun_work(int *, int *, int *, int *);
const int *ocp_model_dyn_impl_ode_fun_sparsity_in(int);
const int *ocp_model_dyn_impl_ode_fun_sparsity_out(int);
int ocp_model_dyn_impl_ode_fun_n_in();
int ocp_model_dyn_impl_ode_fun_n_out();

// implicit ode function jacobian_x_xdot
int ocp_model_dyn_impl_ode_fun_jac_x_xdot(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_dyn_impl_ode_fun_jac_x_xdot_work(int *, int *, int *, int *);
const int *ocp_model_dyn_impl_ode_fun_jac_x_xdot_sparsity_in(int);
const int *ocp_model_dyn_impl_ode_fun_jac_x_xdot_sparsity_out(int);
int ocp_model_dyn_impl_ode_fun_jac_x_xdot_n_in();
int ocp_model_dyn_impl_ode_fun_jac_x_xdot_n_out();

// implicit ode jacobian_x_xdot_u
int ocp_model_dyn_impl_ode_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_dyn_impl_ode_jac_x_xdot_u_work(int *, int *, int *, int *);
const int *ocp_model_dyn_impl_ode_jac_x_xdot_u_sparsity_in(int);
const int *ocp_model_dyn_impl_ode_jac_x_xdot_u_sparsity_out(int);
int ocp_model_dyn_impl_ode_jac_x_xdot_u_n_in();
int ocp_model_dyn_impl_ode_jac_x_xdot_u_n_out();

// implicit ode hessian
int ocp_model_dyn_impl_ode_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_dyn_impl_ode_hess_work(int *, int *, int *, int *);
const int *ocp_model_dyn_impl_ode_hess_sparsity_in(int);
const int *ocp_model_dyn_impl_ode_hess_sparsity_out(int);
int ocp_model_dyn_impl_ode_hess_n_in();
int ocp_model_dyn_impl_ode_hess_n_out();

// discrete model phi function jacobian
int ocp_model_dyn_disc_phi_fun_jac(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_dyn_disc_phi_fun_jac_work(int *, int *, int *, int *);
const int *ocp_model_dyn_disc_phi_fun_jac_sparsity_in(int);
const int *ocp_model_dyn_disc_phi_fun_jac_sparsity_out(int);
int ocp_model_dyn_disc_phi_fun_jac_n_in();
int ocp_model_dyn_disc_phi_fun_jac_n_out();

// discrete model phi function jacobian hessian
int ocp_model_dyn_disc_phi_fun_jac_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_dyn_disc_phi_fun_jac_hess_work(int *, int *, int *, int *);
const int *ocp_model_dyn_disc_phi_fun_jac_hess_sparsity_in(int);
const int *ocp_model_dyn_disc_phi_fun_jac_hess_sparsity_out(int);
int ocp_model_dyn_disc_phi_fun_jac_hess_n_in();
int ocp_model_dyn_disc_phi_fun_jac_hess_n_out();

// nonlinear constraints h
int ocp_model_constr_h_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_constr_h_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *ocp_model_constr_h_fun_jac_ut_xt_sparsity_in(int);
const int *ocp_model_constr_h_fun_jac_ut_xt_sparsity_out(int);
int ocp_model_constr_h_fun_jac_ut_xt_n_in();
int ocp_model_constr_h_fun_jac_ut_xt_n_out();

int ocp_model_constr_h_fun_jac_ut_xt_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_constr_h_fun_jac_ut_xt_hess_work(int *, int *, int *, int *);
const int *ocp_model_constr_h_fun_jac_ut_xt_hess_sparsity_in(int);
const int *ocp_model_constr_h_fun_jac_ut_xt_hess_sparsity_out(int);
int ocp_model_constr_h_fun_jac_ut_xt_hess_n_in();
int ocp_model_constr_h_fun_jac_ut_xt_hess_n_out();

// nonlinear constraints h_e
int ocp_model_constr_h_e_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_constr_h_e_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *ocp_model_constr_h_e_fun_jac_ut_xt_sparsity_in(int);
const int *ocp_model_constr_h_e_fun_jac_ut_xt_sparsity_out(int);
int ocp_model_constr_h_e_fun_jac_ut_xt_n_in();
int ocp_model_constr_h_e_fun_jac_ut_xt_n_out();

int ocp_model_constr_h_e_fun_jac_ut_xt_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_constr_h_e_fun_jac_ut_xt_hess_work(int *, int *, int *, int *);
const int *ocp_model_constr_h_e_fun_jac_ut_xt_hess_sparsity_in(int);
const int *ocp_model_constr_h_e_fun_jac_ut_xt_hess_sparsity_out(int);
int ocp_model_constr_h_e_fun_jac_ut_xt_hess_n_in();
int ocp_model_constr_h_e_fun_jac_ut_xt_hess_n_out();

// nonlinear least squares y
int ocp_model_cost_y_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_cost_y_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *ocp_model_cost_y_fun_jac_ut_xt_sparsity_in(int);
const int *ocp_model_cost_y_fun_jac_ut_xt_sparsity_out(int);
int ocp_model_cost_y_fun_jac_ut_xt_n_in();
int ocp_model_cost_y_fun_jac_ut_xt_n_out();

int ocp_model_cost_y_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_cost_y_hess_work(int *, int *, int *, int *);
const int *ocp_model_cost_y_hess_sparsity_in(int);
const int *ocp_model_cost_y_hess_sparsity_out(int);
int ocp_model_cost_y_hess_n_in();
int ocp_model_cost_y_hess_n_out();

// nonlinear least squares y_e
int ocp_model_cost_y_e_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_cost_y_e_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *ocp_model_cost_y_e_fun_jac_ut_xt_sparsity_in(int);
const int *ocp_model_cost_y_e_fun_jac_ut_xt_sparsity_out(int);
int ocp_model_cost_y_e_fun_jac_ut_xt_n_in();
int ocp_model_cost_y_e_fun_jac_ut_xt_n_out();

int ocp_model_cost_y_e_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_cost_y_e_hess_work(int *, int *, int *, int *);
const int *ocp_model_cost_y_e_hess_sparsity_in(int);
const int *ocp_model_cost_y_e_hess_sparsity_out(int);
int ocp_model_cost_y_e_hess_n_in();
int ocp_model_cost_y_e_hess_n_out();

// external cost
int ocp_model_cost_ext_cost_jac_hes(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_cost_ext_cost_jac_hes_work(int *, int *, int *, int *);
const int *ocp_model_cost_ext_cost_jac_hes_sparsity_in(int);
const int *ocp_model_cost_ext_cost_jac_hes_sparsity_out(int);
int ocp_model_cost_ext_cost_jac_hes_n_in();
int ocp_model_cost_ext_cost_jac_hes_n_out();

// external cost e
int ocp_model_cost_ext_cost_e_jac_hes(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_cost_ext_cost_e_jac_hes_work(int *, int *, int *, int *);
const int *ocp_model_cost_ext_cost_e_jac_hes_sparsity_in(int);
const int *ocp_model_cost_ext_cost_e_jac_hes_sparsity_out(int);
int ocp_model_cost_ext_cost_e_jac_hes_n_in();
int ocp_model_cost_ext_cost_e_jac_hes_n_out();

