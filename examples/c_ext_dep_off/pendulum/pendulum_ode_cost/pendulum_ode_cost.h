/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef pendulum_ode_COST
#define pendulum_ode_COST

#ifdef __cplusplus
extern "C" {
#endif


// Cost at initial shooting node

int pendulum_ode_cost_y_0_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_cost_y_0_fun_work(int *, int *, int *, int *);
const int *pendulum_ode_cost_y_0_fun_sparsity_in(int);
const int *pendulum_ode_cost_y_0_fun_sparsity_out(int);
int pendulum_ode_cost_y_0_fun_n_in(void);
int pendulum_ode_cost_y_0_fun_n_out(void);

int pendulum_ode_cost_y_0_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_cost_y_0_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *pendulum_ode_cost_y_0_fun_jac_ut_xt_sparsity_in(int);
const int *pendulum_ode_cost_y_0_fun_jac_ut_xt_sparsity_out(int);
int pendulum_ode_cost_y_0_fun_jac_ut_xt_n_in(void);
int pendulum_ode_cost_y_0_fun_jac_ut_xt_n_out(void);

int pendulum_ode_cost_y_0_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_cost_y_0_hess_work(int *, int *, int *, int *);
const int *pendulum_ode_cost_y_0_hess_sparsity_in(int);
const int *pendulum_ode_cost_y_0_hess_sparsity_out(int);
int pendulum_ode_cost_y_0_hess_n_in(void);
int pendulum_ode_cost_y_0_hess_n_out(void);



// Cost at path shooting node

int pendulum_ode_cost_y_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_cost_y_fun_work(int *, int *, int *, int *);
const int *pendulum_ode_cost_y_fun_sparsity_in(int);
const int *pendulum_ode_cost_y_fun_sparsity_out(int);
int pendulum_ode_cost_y_fun_n_in(void);
int pendulum_ode_cost_y_fun_n_out(void);

int pendulum_ode_cost_y_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_cost_y_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *pendulum_ode_cost_y_fun_jac_ut_xt_sparsity_in(int);
const int *pendulum_ode_cost_y_fun_jac_ut_xt_sparsity_out(int);
int pendulum_ode_cost_y_fun_jac_ut_xt_n_in(void);
int pendulum_ode_cost_y_fun_jac_ut_xt_n_out(void);

int pendulum_ode_cost_y_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_cost_y_hess_work(int *, int *, int *, int *);
const int *pendulum_ode_cost_y_hess_sparsity_in(int);
const int *pendulum_ode_cost_y_hess_sparsity_out(int);
int pendulum_ode_cost_y_hess_n_in(void);
int pendulum_ode_cost_y_hess_n_out(void);



// Cost at terminal shooting node

int pendulum_ode_cost_y_e_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_cost_y_e_fun_work(int *, int *, int *, int *);
const int *pendulum_ode_cost_y_e_fun_sparsity_in(int);
const int *pendulum_ode_cost_y_e_fun_sparsity_out(int);
int pendulum_ode_cost_y_e_fun_n_in(void);
int pendulum_ode_cost_y_e_fun_n_out(void);

int pendulum_ode_cost_y_e_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_cost_y_e_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *pendulum_ode_cost_y_e_fun_jac_ut_xt_sparsity_in(int);
const int *pendulum_ode_cost_y_e_fun_jac_ut_xt_sparsity_out(int);
int pendulum_ode_cost_y_e_fun_jac_ut_xt_n_in(void);
int pendulum_ode_cost_y_e_fun_jac_ut_xt_n_out(void);

int pendulum_ode_cost_y_e_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_ode_cost_y_e_hess_work(int *, int *, int *, int *);
const int *pendulum_ode_cost_y_e_hess_sparsity_in(int);
const int *pendulum_ode_cost_y_e_hess_sparsity_out(int);
int pendulum_ode_cost_y_e_hess_n_in(void);
int pendulum_ode_cost_y_e_hess_n_out(void);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // pendulum_ode_COST
