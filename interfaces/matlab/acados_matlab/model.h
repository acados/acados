// TODO templetize with casadi function names !!!

// explicit ode function
int model_expl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int model_expl_ode_fun_work(int *, int *, int *, int *);
const int *model_expl_ode_fun_sparsity_in(int);
const int *model_expl_ode_fun_sparsity_out(int);
int model_expl_ode_fun_n_in();
int model_expl_ode_fun_n_out();

// explicit vde forward
int model_expl_vde_for(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int model_expl_vde_for_work(int *, int *, int *, int *);
const int *model_expl_vde_for_sparsity_in(int);
const int *model_expl_vde_for_sparsity_out(int);
int model_expl_vde_for_n_in();
int model_expl_vde_for_n_out();

// implicit ode function
int model_impl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int model_impl_ode_fun_work(int *, int *, int *, int *);
const int *model_impl_ode_fun_sparsity_in(int);
const int *model_impl_ode_fun_sparsity_out(int);
int model_impl_ode_fun_n_in();
int model_impl_ode_fun_n_out();

// implicit ode function jacobian_x_xdot
int model_impl_ode_fun_jac_x_xdot(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int model_impl_ode_fun_jac_x_xdot_work(int *, int *, int *, int *);
const int *model_impl_ode_fun_jac_x_xdot_sparsity_in(int);
const int *model_impl_ode_fun_jac_x_xdot_sparsity_out(int);
int model_impl_ode_fun_jac_x_xdot_n_in();
int model_impl_ode_fun_jac_x_xdot_n_out();

// implicit ode jacobian_x_xdot_u
int model_impl_ode_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int model_impl_ode_jac_x_xdot_u_work(int *, int *, int *, int *);
const int *model_impl_ode_jac_x_xdot_u_sparsity_in(int);
const int *model_impl_ode_jac_x_xdot_u_sparsity_out(int);
int model_impl_ode_jac_x_xdot_u_n_in();
int model_impl_ode_jac_x_xdot_u_n_out();

