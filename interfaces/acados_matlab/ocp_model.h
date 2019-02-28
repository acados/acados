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

// nonlinear constraints h
int ocp_model_constr_h_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_constr_h_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *ocp_model_constr_h_fun_jac_ut_xt_sparsity_in(int);
const int *ocp_model_constr_h_fun_jac_ut_xt_sparsity_out(int);
int ocp_model_constr_h_fun_jac_ut_xt_n_in();
int ocp_model_constr_h_fun_jac_ut_xt_n_out();

// nonlinear constraints h_e
int ocp_model_constr_h_e_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_constr_h_e_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *ocp_model_constr_h_e_fun_jac_ut_xt_sparsity_in(int);
const int *ocp_model_constr_h_e_fun_jac_ut_xt_sparsity_out(int);
int ocp_model_constr_h_e_fun_jac_ut_xt_n_in();
int ocp_model_constr_h_e_fun_jac_ut_xt_n_out();

// nonlinear least squares y
int ocp_model_cost_y_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_cost_y_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *ocp_model_cost_y_fun_jac_ut_xt_sparsity_in(int);
const int *ocp_model_cost_y_fun_jac_ut_xt_sparsity_out(int);
int ocp_model_cost_y_fun_jac_ut_xt_n_in();
int ocp_model_cost_y_fun_jac_ut_xt_n_out();

// nonlinear least squares y_e
int ocp_model_cost_y_e_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int ocp_model_cost_y_e_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *ocp_model_cost_y_e_fun_jac_ut_xt_sparsity_in(int);
const int *ocp_model_cost_y_e_fun_jac_ut_xt_sparsity_out(int);
int ocp_model_cost_y_e_fun_jac_ut_xt_n_in();
int ocp_model_cost_y_e_fun_jac_ut_xt_n_out();

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

