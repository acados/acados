/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
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
