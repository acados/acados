#ifndef EXAMPLES_C_WT_MODEL_NX3_H_
#define EXAMPLES_C_WT_MODEL_NX3_H_

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

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_WT_MODEL_NX3_H_

