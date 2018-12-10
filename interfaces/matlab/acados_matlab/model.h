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

