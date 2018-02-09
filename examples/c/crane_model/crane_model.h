#ifndef EXAMPLES_C_CRANE_MODEL_CRANE_MODEL_H_
#define EXAMPLES_C_CRANE_MODEL_CRANE_MODEL_H_

#include "acados/utils/types.h"

#ifdef __cplusplus
extern "C" {
#endif

int vdeFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
const int* vdeFun_sparsity_out(int i);
int vdeFun_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);

int jacFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
const int* jacFun_sparsity_out(int i);
int jacFun_work(int* sz_arg, int* sz_res, int* sz_iw, int* sz_w);

int adjFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
const int* adjFun_sparsity_out(int i);
int adjFun_work(int* sz_arg, int* sz_res, int* sz_iw, int* sz_w);

int hessFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
const int* hessFun_sparsity_out(int i);
int hessFun_work(int* sz_arg, int* sz_res, int* sz_iw, int* sz_w);

int impl_odeFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
const int* impl_odeFun_sparsity_out(int i);
int impl_odeFun_work(int* sz_arg, int* sz_res, int* sz_iw, int* sz_w);

int impl_jacFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
const int* impl_jacFun_sparsity_out(int i);
int impl_jacFun_work(int* sz_arg, int* sz_res, int* sz_iw, int* sz_w);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_CRANE_MODEL_CRANE_MODEL_H_
