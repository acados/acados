#ifndef EXAMPLES_C_GNSF_CRANE_MODEL_GNSF_CRANE_MODEL_H_
#define EXAMPLES_C_GNSF_CRANE_MODEL_GNSF_CRANE_MODEL_H_

#include "acados/utils/types.h"

#ifdef __cplusplus
extern "C" {
#endif

// used to import integers & double matrices
int get_ints_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int But_KK_YY_ZZ_LO_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

// Phi_inc_dy_fun
int Phi_inc_dy_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int Phi_inc_dy_fun_work(int *, int *, int *, int *);
const int *Phi_inc_dy_fun_sparsity_in(int);
const int *Phi_inc_dy_fun_sparsity_out(int);

// f_LO_inc_J_x1k1uz_fun
int f_LO_inc_J_x1k1uz_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int f_LO_inc_J_x1k1uz_fun_work(int *, int *, int *, int *);
const int *f_LO_inc_J_x1k1uz_fun_sparsity_in(int);
const int *f_LO_inc_J_x1k1uz_fun_sparsity_out(int);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_CRANE_MODEL_CRANE_MODEL_H_