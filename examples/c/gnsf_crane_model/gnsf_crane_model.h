#ifndef EXAMPLES_C_GNSF_CRANE_MODEL_GNSF_CRANE_MODEL_H_
#define EXAMPLES_C_GNSF_CRANE_MODEL_GNSF_CRANE_MODEL_H_

#include "acados/utils/types.h"

#ifdef __cplusplus
extern "C" {
#endif

// int vdeFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int get_ints_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int KKmat_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int ALO_M2_dK2dx2_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int ZZmat_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int Butcher_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int res_inc_Jff_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int jac_res_ffx1u_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int f_LO_inc_J_x1k1uz_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_CRANE_MODEL_CRANE_MODEL_H_