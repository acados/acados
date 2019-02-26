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

