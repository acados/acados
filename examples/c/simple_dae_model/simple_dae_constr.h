/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef EXAMPLES_C_SIMPLE_DAE_CONSTR
#define EXAMPLES_C_SIMPLE_DAE_CONSTR

#ifdef __cplusplus
extern "C" {
#endif

int simple_dae_constr_h_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int simple_dae_constr_h_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *simple_dae_constr_h_fun_jac_ut_xt_sparsity_in(int);
const int *simple_dae_constr_h_fun_jac_ut_xt_sparsity_out(int);
int simple_dae_constr_h_fun_jac_ut_xt_n_in();
int simple_dae_constr_h_fun_jac_ut_xt_n_out();

int        simple_dae_constr_h_fun_jac_ut_xt_hess_(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int        simple_dae_constr_h_fun_jac_ut_xt_hess_work(int *, int *, int *, int *);
const int *simple_dae_constr_h_fun_jac_ut_xt_hess_sparsity_in(int);
const int *simple_dae_constr_h_fun_jac_ut_xt_hess_sparsity_out(int);
int        simple_dae_constr_h_fun_jac_ut_xt_hess_n_in();
int        simple_dae_constr_h_fun_jac_ut_xt_hess_n_out();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_SIMPLE_DAE_CONSTR
