#ifndef EXAMPLES_C_CHAIN_MODEL_CHAIN_MODEL_IMPL_H_
#define EXAMPLES_C_CHAIN_MODEL_CHAIN_MODEL_IMPL_H_

#include "acados/utils/types.h"

#ifdef __cplusplus
extern "C" {
#endif

// implicit ODE
int casadi_impl_ode_fun_chain_nm2(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_chain_nm3(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_chain_nm4(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_chain_nm5(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_chain_nm6(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_chain_nm7(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_chain_nm8(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_chain_nm9(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);

int casadi_impl_ode_fun_chain_nm2_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_chain_nm3_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_chain_nm4_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_chain_nm5_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_chain_nm6_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_chain_nm7_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_chain_nm8_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_chain_nm9_work(int *, int *, int *, int *);

const int *casadi_impl_ode_fun_chain_nm2_sparsity_in(int);
const int *casadi_impl_ode_fun_chain_nm3_sparsity_in(int);
const int *casadi_impl_ode_fun_chain_nm4_sparsity_in(int);
const int *casadi_impl_ode_fun_chain_nm5_sparsity_in(int);
const int *casadi_impl_ode_fun_chain_nm6_sparsity_in(int);
const int *casadi_impl_ode_fun_chain_nm7_sparsity_in(int);
const int *casadi_impl_ode_fun_chain_nm8_sparsity_in(int);
const int *casadi_impl_ode_fun_chain_nm9_sparsity_in(int);

const int *casadi_impl_ode_fun_chain_nm2_sparsity_out(int);
const int *casadi_impl_ode_fun_chain_nm3_sparsity_out(int);
const int *casadi_impl_ode_fun_chain_nm4_sparsity_out(int);
const int *casadi_impl_ode_fun_chain_nm5_sparsity_out(int);
const int *casadi_impl_ode_fun_chain_nm6_sparsity_out(int);
const int *casadi_impl_ode_fun_chain_nm7_sparsity_out(int);
const int *casadi_impl_ode_fun_chain_nm8_sparsity_out(int);
const int *casadi_impl_ode_fun_chain_nm9_sparsity_out(int);

int casadi_impl_ode_fun_chain_nm2_n_in();
int casadi_impl_ode_fun_chain_nm3_n_in();
int casadi_impl_ode_fun_chain_nm4_n_in();
int casadi_impl_ode_fun_chain_nm5_n_in();
int casadi_impl_ode_fun_chain_nm6_n_in();
int casadi_impl_ode_fun_chain_nm7_n_in();
int casadi_impl_ode_fun_chain_nm8_n_in();
int casadi_impl_ode_fun_chain_nm9_n_in();

int casadi_impl_ode_fun_chain_nm2_n_out();
int casadi_impl_ode_fun_chain_nm3_n_out();
int casadi_impl_ode_fun_chain_nm4_n_out();
int casadi_impl_ode_fun_chain_nm5_n_out();
int casadi_impl_ode_fun_chain_nm6_n_out();
int casadi_impl_ode_fun_chain_nm7_n_out();
int casadi_impl_ode_fun_chain_nm8_n_out();
int casadi_impl_ode_fun_chain_nm9_n_out();



// impl_ode_fun_jac_x_xdot
int casadi_impl_ode_fun_jac_x_xdot_chain_nm2(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_chain_nm3(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_chain_nm4(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_chain_nm5(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_chain_nm6(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_chain_nm7(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_chain_nm8(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_chain_nm9(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);

int casadi_impl_ode_fun_jac_x_xdot_chain_nm2_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_jac_x_xdot_chain_nm3_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_jac_x_xdot_chain_nm4_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_jac_x_xdot_chain_nm5_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_jac_x_xdot_chain_nm6_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_jac_x_xdot_chain_nm7_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_jac_x_xdot_chain_nm8_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_jac_x_xdot_chain_nm9_work(int *, int *, int *, int *);

const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm2_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm3_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm4_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm5_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm6_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm7_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm8_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm9_sparsity_in(int);

const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm2_sparsity_out(int);
const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm3_sparsity_out(int);
const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm4_sparsity_out(int);
const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm5_sparsity_out(int);
const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm6_sparsity_out(int);
const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm7_sparsity_out(int);
const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm8_sparsity_out(int);
const int *casadi_impl_ode_fun_jac_x_xdot_chain_nm9_sparsity_out(int);

int casadi_impl_ode_fun_jac_x_xdot_chain_nm2_n_in();
int casadi_impl_ode_fun_jac_x_xdot_chain_nm3_n_in();
int casadi_impl_ode_fun_jac_x_xdot_chain_nm4_n_in();
int casadi_impl_ode_fun_jac_x_xdot_chain_nm5_n_in();
int casadi_impl_ode_fun_jac_x_xdot_chain_nm6_n_in();
int casadi_impl_ode_fun_jac_x_xdot_chain_nm7_n_in();
int casadi_impl_ode_fun_jac_x_xdot_chain_nm8_n_in();
int casadi_impl_ode_fun_jac_x_xdot_chain_nm9_n_in();

int casadi_impl_ode_fun_jac_x_xdot_chain_nm2_n_out();
int casadi_impl_ode_fun_jac_x_xdot_chain_nm3_n_out();
int casadi_impl_ode_fun_jac_x_xdot_chain_nm4_n_out();
int casadi_impl_ode_fun_jac_x_xdot_chain_nm5_n_out();
int casadi_impl_ode_fun_jac_x_xdot_chain_nm6_n_out();
int casadi_impl_ode_fun_jac_x_xdot_chain_nm7_n_out();
int casadi_impl_ode_fun_jac_x_xdot_chain_nm8_n_out();
int casadi_impl_ode_fun_jac_x_xdot_chain_nm9_n_out();

//
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm7(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm8(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm9(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);

int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm7_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm8_work(int *, int *, int *, int *);
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm9_work(int *, int *, int *, int *);

const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm7_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm8_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm9_sparsity_in(int);

const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2_sparsity_out(int);
const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3_sparsity_out(int);
const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4_sparsity_out(int);
const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5_sparsity_out(int);
const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_sparsity_out(int);
const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm7_sparsity_out(int);
const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm8_sparsity_out(int);
const int *casadi_impl_ode_fun_jac_x_xdot_u_chain_nm9_sparsity_out(int);

int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2_n_in();
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3_n_in();
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4_n_in();
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5_n_in();
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_n_in();
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm7_n_in();
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm8_n_in();
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm9_n_in();

int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2_n_out();
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3_n_out();
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4_n_out();
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5_n_out();
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_n_out();
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm7_n_out();
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm8_n_out();
int casadi_impl_ode_fun_jac_x_xdot_u_chain_nm9_n_out();

//
int casadi_impl_ode_jac_x_xdot_u_chain_nm2(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_jac_x_xdot_u_chain_nm3(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_jac_x_xdot_u_chain_nm4(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_jac_x_xdot_u_chain_nm5(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_jac_x_xdot_u_chain_nm6(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_jac_x_xdot_u_chain_nm7(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_jac_x_xdot_u_chain_nm8(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_jac_x_xdot_u_chain_nm9(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);

int casadi_impl_ode_jac_x_xdot_u_chain_nm2_work(int *, int *, int *, int *);
int casadi_impl_ode_jac_x_xdot_u_chain_nm3_work(int *, int *, int *, int *);
int casadi_impl_ode_jac_x_xdot_u_chain_nm4_work(int *, int *, int *, int *);
int casadi_impl_ode_jac_x_xdot_u_chain_nm5_work(int *, int *, int *, int *);
int casadi_impl_ode_jac_x_xdot_u_chain_nm6_work(int *, int *, int *, int *);
int casadi_impl_ode_jac_x_xdot_u_chain_nm7_work(int *, int *, int *, int *);
int casadi_impl_ode_jac_x_xdot_u_chain_nm8_work(int *, int *, int *, int *);
int casadi_impl_ode_jac_x_xdot_u_chain_nm9_work(int *, int *, int *, int *);

const int *casadi_impl_ode_jac_x_xdot_u_chain_nm2_sparsity_in(int);
const int *casadi_impl_ode_jac_x_xdot_u_chain_nm3_sparsity_in(int);
const int *casadi_impl_ode_jac_x_xdot_u_chain_nm4_sparsity_in(int);
const int *casadi_impl_ode_jac_x_xdot_u_chain_nm5_sparsity_in(int);
const int *casadi_impl_ode_jac_x_xdot_u_chain_nm6_sparsity_in(int);
const int *casadi_impl_ode_jac_x_xdot_u_chain_nm7_sparsity_in(int);
const int *casadi_impl_ode_jac_x_xdot_u_chain_nm8_sparsity_in(int);
const int *casadi_impl_ode_jac_x_xdot_u_chain_nm9_sparsity_in(int);

const int *casadi_impl_ode_jac_x_xdot_u_chain_nm2_sparsity_out(int);
const int *casadi_impl_ode_jac_x_xdot_u_chain_nm3_sparsity_out(int);
const int *casadi_impl_ode_jac_x_xdot_u_chain_nm4_sparsity_out(int);
const int *casadi_impl_ode_jac_x_xdot_u_chain_nm5_sparsity_out(int);
const int *casadi_impl_ode_jac_x_xdot_u_chain_nm6_sparsity_out(int);
const int *casadi_impl_ode_jac_x_xdot_u_chain_nm7_sparsity_out(int);
const int *casadi_impl_ode_jac_x_xdot_u_chain_nm8_sparsity_out(int);
const int *casadi_impl_ode_jac_x_xdot_u_chain_nm9_sparsity_out(int);

int casadi_impl_ode_jac_x_xdot_u_chain_nm2_n_in();
int casadi_impl_ode_jac_x_xdot_u_chain_nm3_n_in();
int casadi_impl_ode_jac_x_xdot_u_chain_nm4_n_in();
int casadi_impl_ode_jac_x_xdot_u_chain_nm5_n_in();
int casadi_impl_ode_jac_x_xdot_u_chain_nm6_n_in();
int casadi_impl_ode_jac_x_xdot_u_chain_nm7_n_in();
int casadi_impl_ode_jac_x_xdot_u_chain_nm8_n_in();
int casadi_impl_ode_jac_x_xdot_u_chain_nm9_n_in();

int casadi_impl_ode_jac_x_xdot_u_chain_nm2_n_out();
int casadi_impl_ode_jac_x_xdot_u_chain_nm3_n_out();
int casadi_impl_ode_jac_x_xdot_u_chain_nm4_n_out();
int casadi_impl_ode_jac_x_xdot_u_chain_nm5_n_out();
int casadi_impl_ode_jac_x_xdot_u_chain_nm6_n_out();
int casadi_impl_ode_jac_x_xdot_u_chain_nm7_n_out();
int casadi_impl_ode_jac_x_xdot_u_chain_nm8_n_out();
int casadi_impl_ode_jac_x_xdot_u_chain_nm9_n_out();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_CHAIN_MODEL_CHAIN_MODEL_IMPL_H_
