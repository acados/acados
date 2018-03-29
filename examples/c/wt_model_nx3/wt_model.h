#ifndef EXAMPLES_C_WT_MODEL_NX3_H_
#define EXAMPLES_C_WT_MODEL_NX3_H_

#include "acados/utils/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* explicit ODE */

// explicit ODE
int ode_energy_balanced_model(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int ode_energy_balanced_model_work(int *, int *, int *, int *);
const int *ode_energy_balanced_model_sparsity_in(int);
const int *ode_energy_balanced_model_sparsity_out(int);
int ode_energy_balanced_model_n_in();
int ode_energy_balanced_model_n_out();
// forward explicit VDE
int vde_energy_balanced_model(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int vde_energy_balanced_model_work(int *, int *, int *, int *);
const int *vde_energy_balanced_model_sparsity_in(int);
const int *vde_energy_balanced_model_sparsity_out(int);
int vde_energy_balanced_model_n_in();
int vde_energy_balanced_model_n_out();
// adjoing explicit VDE
int vde_adj_energy_balanced_model(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int vde_adj_energy_balanced_model_work(int *, int *, int *, int *);
const int *vde_adj_energy_balanced_model_sparsity_in(int);
const int *vde_adj_energy_balanced_model_sparsity_out(int);
int vde_adj_energy_balanced_model_n_in();
int vde_adj_energy_balanced_model_n_out();
// hessian explicit ODE
int vde_hess_energy_balanced_model(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int vde_hess_energy_balanced_model_work(int *, int *, int *, int *);
const int *vde_hess_energy_balanced_model_sparsity_in(int);
const int *vde_hess_energy_balanced_model_sparsity_out(int);
int vde_hess_energy_balanced_model_n_in();
int vde_hess_energy_balanced_model_n_out();
// hessian explicit ODE
int jac_energy_balanced_model(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int jac_energy_balanced_model_work(int *, int *, int *, int *);
const int *jac_energy_balanced_model_sparsity_in(int);
const int *jac_energy_balanced_model_sparsity_out(int);
int jac_energy_balanced_model_n_in();
int jac_energy_balanced_model_n_out();

/* implicit ODE */

// implicit ODE
int impl_odeFun_energy_balanced_model(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_odeFun_energy_balanced_model_work(int *, int *, int *, int *);
const int *impl_odeFun_energy_balanced_model_sparsity_in(int);
const int *impl_odeFun_energy_balanced_model_sparsity_out(int);
int impl_odeFun_energy_balanced_model_n_in();
int impl_odeFun_energy_balanced_model_n_out();

// jac_x implicit ODE
int impl_jacFun_x_energy_balanced_model(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_x_energy_balanced_model_work(int *, int *, int *, int *);
const int *impl_jacFun_x_energy_balanced_model_sparsity_in(int);
const int *impl_jacFun_x_energy_balanced_model_sparsity_out(int);
int impl_jacFun_x_energy_balanced_model_n_in();
int impl_jacFun_x_energy_balanced_model_n_out();

// jax_xdot implicit ODE
int impl_jacFun_xdot_energy_balanced_model(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_xdot_energy_balanced_model_work(int *, int *, int *, int *);
const int *impl_jacFun_xdot_energy_balanced_model_sparsity_in(int);
const int *impl_jacFun_xdot_energy_balanced_model_sparsity_out(int);
int impl_jacFun_xdot_energy_balanced_model_n_in();
int impl_jacFun_xdot_energy_balanced_model_n_out();

// jax_u implicit ODE
int impl_jacFun_u_energy_balanced_model(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_u_energy_balanced_model_work(int *, int *, int *, int *);
const int *impl_jacFun_u_energy_balanced_model_sparsity_in(int);
const int *impl_jacFun_u_energy_balanced_model_sparsity_out(int);
int impl_jacFun_u_energy_balanced_model_n_in();
int impl_jacFun_u_energy_balanced_model_n_out();



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_WT_MODEL_NX3_H_

