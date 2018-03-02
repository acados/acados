#ifndef EXAMPLES_C_CHAIN_MODEL_CHAIN_MODEL_IMPL_H_
#define EXAMPLES_C_CHAIN_MODEL_CHAIN_MODEL_IMPL_H_

#include "acados/utils/types.h"

#ifdef __cplusplus
extern "C" {
#endif

// implicit ODE
int impl_odeFun_chain_nm2(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_odeFun_chain_nm3(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_odeFun_chain_nm4(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_odeFun_chain_nm5(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_odeFun_chain_nm6(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_odeFun_chain_nm7(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_odeFun_chain_nm8(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_odeFun_chain_nm9(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

int impl_odeFun_chain_nm2_work(int *, int *, int *, int *);
int impl_odeFun_chain_nm3_work(int *, int *, int *, int *);
int impl_odeFun_chain_nm4_work(int *, int *, int *, int *);
int impl_odeFun_chain_nm5_work(int *, int *, int *, int *);
int impl_odeFun_chain_nm6_work(int *, int *, int *, int *);
int impl_odeFun_chain_nm7_work(int *, int *, int *, int *);
int impl_odeFun_chain_nm8_work(int *, int *, int *, int *);
int impl_odeFun_chain_nm9_work(int *, int *, int *, int *);

const int *impl_odeFun_chain_nm2_sparsity_in(int);
const int *impl_odeFun_chain_nm3_sparsity_in(int);
const int *impl_odeFun_chain_nm4_sparsity_in(int);
const int *impl_odeFun_chain_nm5_sparsity_in(int);
const int *impl_odeFun_chain_nm6_sparsity_in(int);
const int *impl_odeFun_chain_nm7_sparsity_in(int);
const int *impl_odeFun_chain_nm8_sparsity_in(int);
const int *impl_odeFun_chain_nm9_sparsity_in(int);

const int *impl_odeFun_chain_nm2_sparsity_out(int);
const int *impl_odeFun_chain_nm3_sparsity_out(int);
const int *impl_odeFun_chain_nm4_sparsity_out(int);
const int *impl_odeFun_chain_nm5_sparsity_out(int);
const int *impl_odeFun_chain_nm6_sparsity_out(int);
const int *impl_odeFun_chain_nm7_sparsity_out(int);
const int *impl_odeFun_chain_nm8_sparsity_out(int);
const int *impl_odeFun_chain_nm9_sparsity_out(int);

// jac_x implicit ODE
int impl_jacFun_x_chain_nm2(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_x_chain_nm3(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_x_chain_nm4(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_x_chain_nm5(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_x_chain_nm6(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_x_chain_nm7(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_x_chain_nm8(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_x_chain_nm9(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

int impl_jacFun_x_chain_nm2_work(int *, int *, int *, int *);
int impl_jacFun_x_chain_nm3_work(int *, int *, int *, int *);
int impl_jacFun_x_chain_nm4_work(int *, int *, int *, int *);
int impl_jacFun_x_chain_nm5_work(int *, int *, int *, int *);
int impl_jacFun_x_chain_nm6_work(int *, int *, int *, int *);
int impl_jacFun_x_chain_nm7_work(int *, int *, int *, int *);
int impl_jacFun_x_chain_nm8_work(int *, int *, int *, int *);
int impl_jacFun_x_chain_nm9_work(int *, int *, int *, int *);

const int *impl_jacFun_x_chain_nm2_sparsity_in(int);
const int *impl_jacFun_x_chain_nm3_sparsity_in(int);
const int *impl_jacFun_x_chain_nm4_sparsity_in(int);
const int *impl_jacFun_x_chain_nm5_sparsity_in(int);
const int *impl_jacFun_x_chain_nm6_sparsity_in(int);
const int *impl_jacFun_x_chain_nm7_sparsity_in(int);
const int *impl_jacFun_x_chain_nm8_sparsity_in(int);
const int *impl_jacFun_x_chain_nm9_sparsity_in(int);

const int *impl_jacFun_x_chain_nm2_sparsity_out(int);
const int *impl_jacFun_x_chain_nm3_sparsity_out(int);
const int *impl_jacFun_x_chain_nm4_sparsity_out(int);
const int *impl_jacFun_x_chain_nm5_sparsity_out(int);
const int *impl_jacFun_x_chain_nm6_sparsity_out(int);
const int *impl_jacFun_x_chain_nm7_sparsity_out(int);
const int *impl_jacFun_x_chain_nm8_sparsity_out(int);
const int *impl_jacFun_x_chain_nm9_sparsity_out(int);

// jax_xdot implicit ODE
int impl_jacFun_xdot_chain_nm2(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_xdot_chain_nm3(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_xdot_chain_nm4(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_xdot_chain_nm5(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_xdot_chain_nm6(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_xdot_chain_nm7(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_xdot_chain_nm8(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_xdot_chain_nm9(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

int impl_jacFun_xdot_chain_nm2_work(int *, int *, int *, int *);
int impl_jacFun_xdot_chain_nm3_work(int *, int *, int *, int *);
int impl_jacFun_xdot_chain_nm4_work(int *, int *, int *, int *);
int impl_jacFun_xdot_chain_nm5_work(int *, int *, int *, int *);
int impl_jacFun_xdot_chain_nm6_work(int *, int *, int *, int *);
int impl_jacFun_xdot_chain_nm7_work(int *, int *, int *, int *);
int impl_jacFun_xdot_chain_nm8_work(int *, int *, int *, int *);
int impl_jacFun_xdot_chain_nm9_work(int *, int *, int *, int *);

const int *impl_jacFun_xdot_chain_nm2_sparsity_in(int);
const int *impl_jacFun_xdot_chain_nm3_sparsity_in(int);
const int *impl_jacFun_xdot_chain_nm4_sparsity_in(int);
const int *impl_jacFun_xdot_chain_nm5_sparsity_in(int);
const int *impl_jacFun_xdot_chain_nm6_sparsity_in(int);
const int *impl_jacFun_xdot_chain_nm7_sparsity_in(int);
const int *impl_jacFun_xdot_chain_nm8_sparsity_in(int);
const int *impl_jacFun_xdot_chain_nm9_sparsity_in(int);

const int *impl_jacFun_xdot_chain_nm2_sparsity_out(int);
const int *impl_jacFun_xdot_chain_nm3_sparsity_out(int);
const int *impl_jacFun_xdot_chain_nm4_sparsity_out(int);
const int *impl_jacFun_xdot_chain_nm5_sparsity_out(int);
const int *impl_jacFun_xdot_chain_nm6_sparsity_out(int);
const int *impl_jacFun_xdot_chain_nm7_sparsity_out(int);
const int *impl_jacFun_xdot_chain_nm8_sparsity_out(int);
const int *impl_jacFun_xdot_chain_nm9_sparsity_out(int);

// jax_u implicit ODE
int impl_jacFun_u_chain_nm2(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_u_chain_nm3(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_u_chain_nm4(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_u_chain_nm5(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_u_chain_nm6(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_u_chain_nm7(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_u_chain_nm8(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_u_chain_nm9(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

int impl_jacFun_u_chain_nm2_work(int *, int *, int *, int *);
int impl_jacFun_u_chain_nm3_work(int *, int *, int *, int *);
int impl_jacFun_u_chain_nm4_work(int *, int *, int *, int *);
int impl_jacFun_u_chain_nm5_work(int *, int *, int *, int *);
int impl_jacFun_u_chain_nm6_work(int *, int *, int *, int *);
int impl_jacFun_u_chain_nm7_work(int *, int *, int *, int *);
int impl_jacFun_u_chain_nm8_work(int *, int *, int *, int *);
int impl_jacFun_u_chain_nm9_work(int *, int *, int *, int *);

const int *impl_jacFun_u_chain_nm2_sparsity_in(int);
const int *impl_jacFun_u_chain_nm3_sparsity_in(int);
const int *impl_jacFun_u_chain_nm4_sparsity_in(int);
const int *impl_jacFun_u_chain_nm5_sparsity_in(int);
const int *impl_jacFun_u_chain_nm6_sparsity_in(int);
const int *impl_jacFun_u_chain_nm7_sparsity_in(int);
const int *impl_jacFun_u_chain_nm8_sparsity_in(int);
const int *impl_jacFun_u_chain_nm9_sparsity_in(int);

const int *impl_jacFun_u_chain_nm2_sparsity_out(int);
const int *impl_jacFun_u_chain_nm3_sparsity_out(int);
const int *impl_jacFun_u_chain_nm4_sparsity_out(int);
const int *impl_jacFun_u_chain_nm5_sparsity_out(int);
const int *impl_jacFun_u_chain_nm6_sparsity_out(int);
const int *impl_jacFun_u_chain_nm7_sparsity_out(int);
const int *impl_jacFun_u_chain_nm8_sparsity_out(int);
const int *impl_jacFun_u_chain_nm9_sparsity_out(int);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_CHAIN_MODEL_CHAIN_MODEL_IMPL_H_
