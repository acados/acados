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


#ifndef EXAMPLES_C_INV_PENDULUM
#define EXAMPLES_C_INV_PENDULUM

#ifdef __cplusplus
extern "C" {
#endif

// this is a crane model with an artificially added algebraic equation to test gnsf & dae integrators
/* implicit ODE */

// implicit ODE
int pendulum_dae_dyn_impl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_dae_dyn_impl_ode_fun_work(int *, int *, int *, int *);
const int *pendulum_dae_dyn_impl_ode_fun_sparsity_in(int);
const int *pendulum_dae_dyn_impl_ode_fun_sparsity_out(int);
int pendulum_dae_dyn_impl_ode_fun_n_in();
int pendulum_dae_dyn_impl_ode_fun_n_out();

// implicit ODE
int pendulum_dae_dyn_impl_ode_fun_jac_x_xdot(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_dae_dyn_impl_ode_fun_jac_x_xdot_work(int *, int *, int *, int *);
const int *pendulum_dae_dyn_impl_ode_fun_jac_x_xdot_sparsity_in(int);
const int *pendulum_dae_dyn_impl_ode_fun_jac_x_xdot_sparsity_out(int);
int pendulum_dae_dyn_impl_ode_fun_jac_x_xdot_n_in();
int pendulum_dae_dyn_impl_ode_fun_jac_x_xdot_n_out();

// implicit ODE
int pendulum_dae_dyn_impl_ode_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_dae_dyn_impl_ode_jac_x_xdot_u_work(int *, int *, int *, int *);
const int *pendulum_dae_dyn_impl_ode_jac_x_xdot_u_sparsity_in(int);
const int *pendulum_dae_dyn_impl_ode_jac_x_xdot_u_sparsity_out(int);
int pendulum_dae_dyn_impl_ode_jac_x_xdot_u_n_in();
int pendulum_dae_dyn_impl_ode_jac_x_xdot_u_n_out();

// implicit ODE - for lifted_irk
int pendulum_dae_dyn_impl_ode_fun_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_dae_dyn_impl_ode_fun_jac_x_xdot_u_work(int *, int *, int *, int *);
const int *pendulum_dae_dyn_impl_ode_fun_jac_x_xdot_u_sparsity_in(int);
const int *pendulum_dae_dyn_impl_ode_fun_jac_x_xdot_u_sparsity_out(int);
int pendulum_dae_dyn_impl_ode_fun_jac_x_xdot_u_n_in();
int pendulum_dae_dyn_impl_ode_fun_jac_x_xdot_u_n_out();

// implicit ODE
int pendulum_dae_dyn_impl_ode_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int pendulum_dae_dyn_impl_ode_hess_work(int *, int *, int *, int *);
const int *pendulum_dae_dyn_impl_ode_hess_sparsity_in(int);
const int *pendulum_dae_dyn_impl_ode_hess_sparsity_out(int);
int pendulum_dae_dyn_impl_ode_hess_n_in();
int pendulum_dae_dyn_impl_ode_hess_n_out();

/* GNSF Functions */
// used to import model matrices
int        pendulum_dae_dyn_gnsf_get_matrices_fun(const double** arg, double** res, int* iw, double* w, void *mem);
int        pendulum_dae_dyn_gnsf_get_matrices_fun_work(int *, int *, int *, int *);
const int *pendulum_dae_dyn_gnsf_get_matrices_fun_sparsity_in(int);
const int *pendulum_dae_dyn_gnsf_get_matrices_fun_sparsity_out(int);
int        pendulum_dae_dyn_gnsf_get_matrices_fun_n_in();
int        pendulum_dae_dyn_gnsf_get_matrices_fun_n_out();

// phi_fun
int        pendulum_dae_dyn_gnsf_phi_fun(const double** arg, double** res, int* iw, double* w, void *mem);
int        pendulum_dae_dyn_gnsf_phi_fun_work(int *, int *, int *, int *);
const int *pendulum_dae_dyn_gnsf_phi_fun_sparsity_in(int);
const int *pendulum_dae_dyn_gnsf_phi_fun_sparsity_out(int);
int        pendulum_dae_dyn_gnsf_phi_fun_n_in();
int        pendulum_dae_dyn_gnsf_phi_fun_n_out();

// phi_fun_jac_y
int        pendulum_dae_dyn_gnsf_phi_fun_jac_y(const double** arg, double** res, int* iw, double* w, void *mem);
int        pendulum_dae_dyn_gnsf_phi_fun_jac_y_work(int *, int *, int *, int *);
const int *pendulum_dae_dyn_gnsf_phi_fun_jac_y_sparsity_in(int);
const int *pendulum_dae_dyn_gnsf_phi_fun_jac_y_sparsity_out(int);
int        pendulum_dae_dyn_gnsf_phi_fun_jac_y_n_in();
int        pendulum_dae_dyn_gnsf_phi_fun_jac_y_n_out();

// phi_jac_y_uhat
int        pendulum_dae_dyn_gnsf_phi_jac_y_uhat(const double** arg, double** res, int* iw, double* w, void *mem);
int        pendulum_dae_dyn_gnsf_phi_jac_y_uhat_work(int *, int *, int *, int *);
const int *pendulum_dae_dyn_gnsf_phi_jac_y_uhat_sparsity_in(int);
const int *pendulum_dae_dyn_gnsf_phi_jac_y_uhat_sparsity_out(int);
int        pendulum_dae_dyn_gnsf_phi_jac_y_uhat_n_in();
int        pendulum_dae_dyn_gnsf_phi_jac_y_uhat_n_out();

// f_lo_fun_jac_x1k1uz
int        pendulum_dae_dyn_gnsf_f_lo_fun_jac_x1k1uz(const double** arg, double** res, int* iw, double* w, void *mem);
int        pendulum_dae_dyn_gnsf_f_lo_fun_jac_x1k1uz_work(int *, int *, int *, int *);
const int *pendulum_dae_dyn_gnsf_f_lo_fun_jac_x1k1uz_sparsity_in(int);
const int *pendulum_dae_dyn_gnsf_f_lo_fun_jac_x1k1uz_sparsity_out(int);
int        pendulum_dae_dyn_gnsf_f_lo_fun_jac_x1k1uz_n_in();
int        pendulum_dae_dyn_gnsf_f_lo_fun_jac_x1k1uz_n_out();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_INV_PENDULUM
