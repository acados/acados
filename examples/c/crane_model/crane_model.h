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

#ifndef EXAMPLES_C_CRANE_MODEL_CRANE_MODEL_H_
#define EXAMPLES_C_CRANE_MODEL_CRANE_MODEL_H_

#include "acados/utils/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* explicit ODE */

// forward explicit VDE
int vdeFun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int vdeFun_work(int *, int *, int *, int *);
const int *vdeFun_sparsity_in(int);
const int *vdeFun_sparsity_out(int);
int vdeFun_n_in();
int vdeFun_n_out();
// adjoing explicit VDE
int adjFun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int adjFun_work(int *, int *, int *, int *);
const int *adjFun_sparsity_in(int);
const int *adjFun_sparsity_out(int);
int adjFun_n_in();
int adjFun_n_out();
// hessian explicit ODE
int hessFun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int hessFun_work(int *, int *, int *, int *);
const int *hessFun_sparsity_in(int);
const int *hessFun_sparsity_out(int);
int hessFun_n_in();
int hessFun_n_out();
// hessian explicit ODE
int jacFun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int jacFun_work(int *, int *, int *, int *);
const int *jacFun_sparsity_in(int);
const int *jacFun_sparsity_out(int);
int jacFun_n_in();
int jacFun_n_out();

/* implicit ODE */

// implicit ODE
int casadi_impl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_work(int *, int *, int *, int *);
const int *casadi_impl_ode_fun_sparsity_in(int);
const int *casadi_impl_ode_fun_sparsity_out(int);
int casadi_impl_ode_fun_n_in();
int casadi_impl_ode_fun_n_out();

//
int casadi_impl_ode_fun_jac_x_xdot(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_fun_jac_x_xdot_work(int *, int *, int *, int *);
const int *casadi_impl_ode_fun_jac_x_xdot_sparsity_in(int);
const int *casadi_impl_ode_fun_jac_x_xdot_sparsity_out(int);
int casadi_impl_ode_fun_jac_x_xdot_n_in();
int casadi_impl_ode_fun_jac_x_xdot_n_out();

//
int casadi_impl_ode_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int casadi_impl_ode_jac_x_xdot_u_work(int *, int *, int *, int *);
const int *casadi_impl_ode_jac_x_xdot_u_sparsity_in(int);
const int *casadi_impl_ode_jac_x_xdot_u_sparsity_out(int);
int casadi_impl_ode_jac_x_xdot_u_n_in();
int casadi_impl_ode_jac_x_xdot_u_n_out();


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_CRANE_MODEL_CRANE_MODEL_H_
