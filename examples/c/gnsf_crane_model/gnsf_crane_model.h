#ifndef EXAMPLES_C_GNSF_CRANE_MODEL_GNSF_CRANE_MODEL_H_
#define EXAMPLES_C_GNSF_CRANE_MODEL_GNSF_CRANE_MODEL_H_

#include "acados/utils/types.h"

#ifdef __cplusplus
extern "C" {
#endif

// used to import integers & double matrices
int get_ints_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int But_KK_ZZ_LO_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

// res_inc_J_ff_fun
int res_inc_Jff_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int res_inc_Jff_fun_work(int *, int *, int *, int *);
const int *res_inc_Jff_fun_sparsity_in(int);
const int *res_inc_Jff_fun_sparsity_out(int);
int *res_inc_Jff_fun_n_in();
int *res_inc_Jff_fun_n_out();

// jac_res_ffx1u_fun
int jac_res_ffx1u_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int jac_res_ffx1u_fun_work(int *, int *, int *, int *);
const int *jac_res_ffx1u_fun_sparsity_in(int);
const int *jac_res_ffx1u_fun_sparsity_out(int);
int *jac_res_ffx1u_fun_n_in();
int *jac_res_ffx1u_fun_n_out();

// f_LO_inc_J_x1k1uz_fun
int f_LO_inc_J_x1k1uz_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int f_LO_inc_J_x1k1uz_fun_work(int *, int *, int *, int *);
const int *f_LO_inc_J_x1k1uz_fun_sparsity_in(int);
const int *f_LO_inc_J_x1k1uz_fun_sparsity_out(int);
int *f_LO_inc_J_x1k1uz_fun_n_in();
int *f_LO_inc_J_x1k1uz_fun_n_out();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_CRANE_MODEL_CRANE_MODEL_H_