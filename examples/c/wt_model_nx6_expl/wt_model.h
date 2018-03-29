#ifndef EXAMPLES_C_WT_MODEL_NX3_H_
#define EXAMPLES_C_WT_MODEL_NX3_H_

#ifdef __cplusplus
extern "C" {
#endif

/* explicit ODE */

// explicit ODE
int expl_ode(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int expl_ode_work(int *, int *, int *, int *);
const int *expl_ode_sparsity_in(int);
const int *expl_ode_sparsity_out(int);
int expl_ode_n_in();
int expl_ode_n_out();

// explicit forward VDE
int expl_forw_vde(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int expl_forw_vde_work(int *, int *, int *, int *);
const int *expl_forw_vde_sparsity_in(int);
const int *expl_forw_vde_sparsity_out(int);
int expl_forw_vde_n_in();
int expl_forw_vde_n_out();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_WT_MODEL_NX3_H_

