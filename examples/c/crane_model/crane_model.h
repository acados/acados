#ifndef EXAMPLES_C_CRANE_MODEL_CRANE_MODEL_H_
#define EXAMPLES_C_CRANE_MODEL_CRANE_MODEL_H_

#include "acados/utils/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* explicit ODE */

// forward explicit VDE
int vdeFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int vdeFun_work(int *, int *, int *, int *);
const int *vdeFun_sparsity_in(int);
const int *vdeFun_sparsity_out(int);
// adjoing explicit VDE
int adjFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int adjFun_work(int *, int *, int *, int *);
const int *adjFun_sparsity_in(int);
const int *adjFun_sparsity_out(int);
// hessian explicit ODE
int hessFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int hessFun_work(int *, int *, int *, int *);
const int *hessFun_sparsity_in(int);
const int *hessFun_sparsity_out(int);



int impl_odeFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_x(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_xdot(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_u(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_CRANE_MODEL_CRANE_MODEL_H_
