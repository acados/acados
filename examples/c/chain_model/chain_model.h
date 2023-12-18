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

#ifndef EXAMPLES_CASADI_CHAIN_CHAIN_MODEL_H_
#define EXAMPLES_CASADI_CHAIN_CHAIN_MODEL_H_

#include "acados/utils/types.h"

#ifdef __cplusplus
extern "C" {
#endif


/* forward vde */
int vde_chain_nm2(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_chain_nm3(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_chain_nm4(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_chain_nm5(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_chain_nm6(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_chain_nm7(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_chain_nm8(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_chain_nm9(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);

int vde_chain_nm2_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int vde_chain_nm3_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int vde_chain_nm4_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int vde_chain_nm5_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int vde_chain_nm6_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int vde_chain_nm7_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int vde_chain_nm8_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int vde_chain_nm9_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);

const int* vde_chain_nm2_sparsity_in(int i);
const int* vde_chain_nm3_sparsity_in(int i);
const int* vde_chain_nm4_sparsity_in(int i);
const int* vde_chain_nm5_sparsity_in(int i);
const int* vde_chain_nm6_sparsity_in(int i);
const int* vde_chain_nm7_sparsity_in(int i);
const int* vde_chain_nm8_sparsity_in(int i);
const int* vde_chain_nm9_sparsity_in(int i);

const int* vde_chain_nm2_sparsity_out(int i);
const int* vde_chain_nm3_sparsity_out(int i);
const int* vde_chain_nm4_sparsity_out(int i);
const int* vde_chain_nm5_sparsity_out(int i);
const int* vde_chain_nm6_sparsity_out(int i);
const int* vde_chain_nm7_sparsity_out(int i);
const int* vde_chain_nm8_sparsity_out(int i);
const int* vde_chain_nm9_sparsity_out(int i);

int vde_chain_nm2_n_in();
int vde_chain_nm3_n_in();
int vde_chain_nm4_n_in();
int vde_chain_nm5_n_in();
int vde_chain_nm6_n_in();
int vde_chain_nm7_n_in();
int vde_chain_nm8_n_in();
int vde_chain_nm9_n_in();

int vde_chain_nm2_n_out();
int vde_chain_nm3_n_out();
int vde_chain_nm4_n_out();
int vde_chain_nm5_n_out();
int vde_chain_nm6_n_out();
int vde_chain_nm7_n_out();
int vde_chain_nm8_n_out();
int vde_chain_nm9_n_out();


/* adjoint vde */
int vde_adj_chain_nm2(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_adj_chain_nm3(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_adj_chain_nm4(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_adj_chain_nm5(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_adj_chain_nm6(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_adj_chain_nm7(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_adj_chain_nm8(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_adj_chain_nm9(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);

/* hessian vde (or ode???) */
int vde_hess_chain_nm2(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_hess_chain_nm3(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_hess_chain_nm4(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_hess_chain_nm5(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_hess_chain_nm6(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_hess_chain_nm7(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_hess_chain_nm8(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int vde_hess_chain_nm9(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);

/* ls cost */
int ls_cost_nm2(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int ls_cost_nm3(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int ls_cost_nm4(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int ls_cost_nm5(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int ls_cost_nm6(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int ls_cost_nm7(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int ls_cost_nm8(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int ls_cost_nm9(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);

int ls_cost_nm2_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int ls_cost_nm3_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int ls_cost_nm4_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int ls_cost_nm5_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int ls_cost_nm6_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int ls_cost_nm7_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int ls_cost_nm8_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int ls_cost_nm9_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);

const int* ls_cost_nm2_sparsity_in(int i);
const int* ls_cost_nm3_sparsity_in(int i);
const int* ls_cost_nm4_sparsity_in(int i);
const int* ls_cost_nm5_sparsity_in(int i);
const int* ls_cost_nm6_sparsity_in(int i);
const int* ls_cost_nm7_sparsity_in(int i);
const int* ls_cost_nm8_sparsity_in(int i);
const int* ls_cost_nm9_sparsity_in(int i);

const int* ls_cost_nm2_sparsity_out(int i);
const int* ls_cost_nm3_sparsity_out(int i);
const int* ls_cost_nm4_sparsity_out(int i);
const int* ls_cost_nm5_sparsity_out(int i);
const int* ls_cost_nm6_sparsity_out(int i);
const int* ls_cost_nm7_sparsity_out(int i);
const int* ls_cost_nm8_sparsity_out(int i);
const int* ls_cost_nm9_sparsity_out(int i);

int ls_cost_nm2_n_in();
int ls_cost_nm3_n_in();
int ls_cost_nm4_n_in();
int ls_cost_nm5_n_in();
int ls_cost_nm6_n_in();
int ls_cost_nm7_n_in();
int ls_cost_nm8_n_in();
int ls_cost_nm9_n_in();

int ls_cost_nm2_n_out();
int ls_cost_nm3_n_out();
int ls_cost_nm4_n_out();
int ls_cost_nm5_n_out();
int ls_cost_nm6_n_out();
int ls_cost_nm7_n_out();
int ls_cost_nm8_n_out();
int ls_cost_nm9_n_out();

/* ls cost N */
int ls_costN_nm2(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int ls_costN_nm3(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int ls_costN_nm4(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int ls_costN_nm5(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int ls_costN_nm6(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int ls_costN_nm7(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int ls_costN_nm8(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int ls_costN_nm9(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);

int ls_costN_nm2_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int ls_costN_nm3_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int ls_costN_nm4_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int ls_costN_nm5_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int ls_costN_nm6_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int ls_costN_nm7_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int ls_costN_nm8_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int ls_costN_nm9_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);

const int* ls_costN_nm2_sparsity_in(int i);
const int* ls_costN_nm3_sparsity_in(int i);
const int* ls_costN_nm4_sparsity_in(int i);
const int* ls_costN_nm5_sparsity_in(int i);
const int* ls_costN_nm6_sparsity_in(int i);
const int* ls_costN_nm7_sparsity_in(int i);
const int* ls_costN_nm8_sparsity_in(int i);
const int* ls_costN_nm9_sparsity_in(int i);

const int* ls_costN_nm2_sparsity_out(int i);
const int* ls_costN_nm3_sparsity_out(int i);
const int* ls_costN_nm4_sparsity_out(int i);
const int* ls_costN_nm5_sparsity_out(int i);
const int* ls_costN_nm6_sparsity_out(int i);
const int* ls_costN_nm7_sparsity_out(int i);
const int* ls_costN_nm8_sparsity_out(int i);
const int* ls_costN_nm9_sparsity_out(int i);

int ls_costN_nm2_n_in();
int ls_costN_nm3_n_in();
int ls_costN_nm4_n_in();
int ls_costN_nm5_n_in();
int ls_costN_nm6_n_in();
int ls_costN_nm7_n_in();
int ls_costN_nm8_n_in();
int ls_costN_nm9_n_in();

int ls_costN_nm2_n_out();
int ls_costN_nm3_n_out();
int ls_costN_nm4_n_out();
int ls_costN_nm5_n_out();
int ls_costN_nm6_n_out();
int ls_costN_nm7_n_out();
int ls_costN_nm8_n_out();
int ls_costN_nm9_n_out();


int pathcon_nm2(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int pathcon_nm3(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int pathcon_nm4(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int pathcon_nm5(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int pathcon_nm6(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int pathcon_nm7(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int pathcon_nm8(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int pathcon_nm9(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);

int pathcon_nm2_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int pathcon_nm3_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int pathcon_nm4_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int pathcon_nm5_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int pathcon_nm6_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int pathcon_nm7_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int pathcon_nm8_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int pathcon_nm9_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);

const int* pathcon_nm2_sparsity_out(int i);
const int* pathcon_nm3_sparsity_out(int i);
const int* pathcon_nm4_sparsity_out(int i);
const int* pathcon_nm5_sparsity_out(int i);
const int* pathcon_nm6_sparsity_out(int i);
const int* pathcon_nm7_sparsity_out(int i);
const int* pathcon_nm8_sparsity_out(int i);
const int* pathcon_nm9_sparsity_out(int i);

int pathconN_nm2(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int pathconN_nm3(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int pathconN_nm4(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int pathconN_nm5(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int pathconN_nm6(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int pathconN_nm7(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int pathconN_nm8(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int pathconN_nm9(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);

int pathconN_nm2_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int pathconN_nm3_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int pathconN_nm4_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int pathconN_nm5_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int pathconN_nm6_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int pathconN_nm7_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int pathconN_nm8_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int pathconN_nm9_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);

const int* pathconN_nm2_sparsity_out(int i);
const int* pathconN_nm3_sparsity_out(int i);
const int* pathconN_nm4_sparsity_out(int i);
const int* pathconN_nm5_sparsity_out(int i);
const int* pathconN_nm6_sparsity_out(int i);
const int* pathconN_nm7_sparsity_out(int i);
const int* pathconN_nm8_sparsity_out(int i);
const int* pathconN_nm9_sparsity_out(int i);

/* casadi ERK integrator */
int casadi_erk4_chain_nm2(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int casadi_erk4_chain_nm3(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int casadi_erk4_chain_nm4(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int casadi_erk4_chain_nm5(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int casadi_erk4_chain_nm6(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int casadi_erk4_chain_nm7(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int casadi_erk4_chain_nm8(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int casadi_erk4_chain_nm9(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);

int casadi_erk4_chain_nm2_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int casadi_erk4_chain_nm3_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int casadi_erk4_chain_nm4_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int casadi_erk4_chain_nm5_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int casadi_erk4_chain_nm6_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int casadi_erk4_chain_nm7_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int casadi_erk4_chain_nm8_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int casadi_erk4_chain_nm9_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);

const int* casadi_erk4_chain_nm2_sparsity_in(int i);
const int* casadi_erk4_chain_nm3_sparsity_in(int i);
const int* casadi_erk4_chain_nm4_sparsity_in(int i);
const int* casadi_erk4_chain_nm5_sparsity_in(int i);
const int* casadi_erk4_chain_nm6_sparsity_in(int i);
const int* casadi_erk4_chain_nm7_sparsity_in(int i);
const int* casadi_erk4_chain_nm8_sparsity_in(int i);
const int* casadi_erk4_chain_nm9_sparsity_in(int i);

const int* casadi_erk4_chain_nm2_sparsity_out(int i);
const int* casadi_erk4_chain_nm3_sparsity_out(int i);
const int* casadi_erk4_chain_nm4_sparsity_out(int i);
const int* casadi_erk4_chain_nm5_sparsity_out(int i);
const int* casadi_erk4_chain_nm6_sparsity_out(int i);
const int* casadi_erk4_chain_nm7_sparsity_out(int i);
const int* casadi_erk4_chain_nm8_sparsity_out(int i);
const int* casadi_erk4_chain_nm9_sparsity_out(int i);

int casadi_erk4_chain_nm2_n_in();
int casadi_erk4_chain_nm3_n_in();
int casadi_erk4_chain_nm4_n_in();
int casadi_erk4_chain_nm5_n_in();
int casadi_erk4_chain_nm6_n_in();
int casadi_erk4_chain_nm7_n_in();
int casadi_erk4_chain_nm8_n_in();
int casadi_erk4_chain_nm9_n_in();

int casadi_erk4_chain_nm2_n_out();
int casadi_erk4_chain_nm3_n_out();
int casadi_erk4_chain_nm4_n_out();
int casadi_erk4_chain_nm5_n_out();
int casadi_erk4_chain_nm6_n_out();
int casadi_erk4_chain_nm7_n_out();
int casadi_erk4_chain_nm8_n_out();
int casadi_erk4_chain_nm9_n_out();

/* external cost */
int chain_nm_2_external_cost(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int chain_nm_3_external_cost(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int chain_nm_4_external_cost(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int chain_nm_5_external_cost(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);
int chain_nm_6_external_cost(const real_t **arg, real_t **res, int *iw, real_t *w, void *mem);

int chain_nm_2_external_cost_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int chain_nm_3_external_cost_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int chain_nm_4_external_cost_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int chain_nm_5_external_cost_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
int chain_nm_6_external_cost_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);

const int* chain_nm_2_external_cost_sparsity_in(int i);
const int* chain_nm_3_external_cost_sparsity_in(int i);
const int* chain_nm_4_external_cost_sparsity_in(int i);
const int* chain_nm_5_external_cost_sparsity_in(int i);
const int* chain_nm_6_external_cost_sparsity_in(int i);

const int* chain_nm_2_external_cost_sparsity_out(int i);
const int* chain_nm_3_external_cost_sparsity_out(int i);
const int* chain_nm_4_external_cost_sparsity_out(int i);
const int* chain_nm_5_external_cost_sparsity_out(int i);
const int* chain_nm_6_external_cost_sparsity_out(int i);

int chain_nm_2_external_cost_n_in();
int chain_nm_3_external_cost_n_in();
int chain_nm_4_external_cost_n_in();
int chain_nm_5_external_cost_n_in();
int chain_nm_6_external_cost_n_in();

int chain_nm_2_external_cost_n_out();
int chain_nm_3_external_cost_n_out();
int chain_nm_4_external_cost_n_out();
int chain_nm_5_external_cost_n_out();
int chain_nm_6_external_cost_n_out();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_CASADI_CHAIN_CHAIN_MODEL_H_
