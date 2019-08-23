/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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

#ifndef {{ ocp.model_name }}_MODEL
#define {{ ocp.model_name }}_MODEL

#ifdef __cplusplus
extern "C" {
#endif

{% if ocp.solver_config.integrator_type != "ERK" %}
// implicit ODE
int {{ ocp.model_name }}_impl_dae_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ocp.model_name }}_impl_dae_fun_work(int *, int *, int *, int *);
const int *{{ ocp.model_name }}_impl_dae_fun_sparsity_in(int);
const int *{{ ocp.model_name }}_impl_dae_fun_sparsity_out(int);
int {{ ocp.model_name }}_impl_dae_fun_n_in();
int {{ ocp.model_name }}_impl_dae_fun_n_out();

// implicit ODE
int {{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_z(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_z_work(int *, int *, int *, int *);
const int *{{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_z_sparsity_in(int);
const int *{{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_z_sparsity_out(int);
int {{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_z_n_in();
int {{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_z_n_out();

// implicit ODE
int {{ ocp.model_name }}_impl_dae_jac_x_xdot_u_z(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ocp.model_name }}_impl_dae_jac_x_xdot_u_z_work(int *, int *, int *, int *);
const int *{{ ocp.model_name }}_impl_dae_jac_x_xdot_u_z_sparsity_in(int);
const int *{{ ocp.model_name }}_impl_dae_jac_x_xdot_u_z_sparsity_out(int);
int {{ ocp.model_name }}_impl_dae_jac_x_xdot_u_z_n_in();
int {{ ocp.model_name }}_impl_dae_jac_x_xdot_u_z_n_out();

// // implicit ODE - for lifted_irk
// int {{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_u(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
// int {{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_u_work(int *, int *, int *, int *);
// const int *{{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_u_sparsity_in(int);
// const int *{{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_u_sparsity_out(int);
// int {{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_u_n_in();
// int {{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_u_n_out();

// // implicit ODE
// int {{ ocp.model_name }}_impl_dae_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
// int {{ ocp.model_name }}_impl_dae_hess_work(int *, int *, int *, int *);
// const int *{{ ocp.model_name }}_impl_dae_hess_sparsity_in(int);
// const int *{{ ocp.model_name }}_impl_dae_hess_sparsity_out(int);
// int {{ ocp.model_name }}_impl_dae_hess_n_in();
// int {{ ocp.model_name }}_impl_dae_hess_n_out();

// /* GNSF Functions */
// // used to import model matrices
// int        {{ ocp.model_name }}_get_matrices_fun(const double** arg, double** res, int* iw, double* w, void *mem);
// int        {{ ocp.model_name }}_get_matrices_fun_work(int *, int *, int *, int *);
// const int *{{ ocp.model_name }}_get_matrices_fun_sparsity_in(int);
// const int *{{ ocp.model_name }}_get_matrices_fun_sparsity_out(int);
// int        {{ ocp.model_name }}_get_matrices_fun_n_in();
// int        {{ ocp.model_name }}_get_matrices_fun_n_out();

// // phi_fun
// int        {{ ocp.model_name }}_phi_fun(const double** arg, double** res, int* iw, double* w, void *mem);
// int        {{ ocp.model_name }}_phi_fun_work(int *, int *, int *, int *);
// const int *{{ ocp.model_name }}_phi_fun_sparsity_in(int);
// const int *{{ ocp.model_name }}_phi_fun_sparsity_out(int);
// int        {{ ocp.model_name }}_phi_fun_n_in();
// int        {{ ocp.model_name }}_phi_fun_n_out();

// // phi_fun_jac_y
// int        {{ ocp.model_name }}_phi_fun_jac_y(const double** arg, double** res, int* iw, double* w, void *mem);
// int        {{ ocp.model_name }}_phi_fun_jac_y_work(int *, int *, int *, int *);
// const int *{{ ocp.model_name }}_phi_fun_jac_y_sparsity_in(int);
// const int *{{ ocp.model_name }}_phi_fun_jac_y_sparsity_out(int);
// int        {{ ocp.model_name }}_phi_fun_jac_y_n_in();
// int        {{ ocp.model_name }}_phi_fun_jac_y_n_out();

// // phi_jac_y_uhat
// int        {{ ocp.model_name }}_phi_jac_y_uhat(const double** arg, double** res, int* iw, double* w, void *mem);
// int        {{ ocp.model_name }}_phi_jac_y_uhat_work(int *, int *, int *, int *);
// const int *{{ ocp.model_name }}_phi_jac_y_uhat_sparsity_in(int);
// const int *{{ ocp.model_name }}_phi_jac_y_uhat_sparsity_out(int);
// int        {{ ocp.model_name }}_phi_jac_y_uhat_n_in();
// int        {{ ocp.model_name }}_phi_jac_y_uhat_n_out();

// // f_lo_fun_jac_x1k1uz
// int        {{ ocp.model_name }}_f_lo_fun_jac_x1k1uz(const double** arg, double** res, int* iw, double* w, void *mem);
// int        {{ ocp.model_name }}_f_lo_fun_jac_x1k1uz_work(int *, int *, int *, int *);
// const int *{{ ocp.model_name }}_f_lo_fun_jac_x1k1uz_sparsity_in(int);
// const int *{{ ocp.model_name }}_f_lo_fun_jac_x1k1uz_sparsity_out(int);
// int        {{ ocp.model_name }}_f_lo_fun_jac_x1k1uz_n_in();
// int        {{ ocp.model_name }}_f_lo_fun_jac_x1k1uz_n_out();

{% else %}
/* explicit ODE */

// explicit ODE
int {{ ocp.model_name }}_expl_ode_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ocp.model_name }}_expl_ode_fun_work(int *, int *, int *, int *);
const int *{{ ocp.model_name }}_expl_ode_fun_sparsity_in(int);
const int *{{ ocp.model_name }}_expl_ode_fun_sparsity_out(int);
int {{ ocp.model_name }}_expl_ode_fun_n_in();
int {{ ocp.model_name }}_expl_ode_fun_n_out();

// explicit forward VDE
int {{ ocp.model_name }}_expl_vde_forw(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ocp.model_name }}_expl_vde_forw_work(int *, int *, int *, int *);
const int *{{ ocp.model_name }}_expl_vde_forw_sparsity_in(int);
const int *{{ ocp.model_name }}_expl_vde_forw_sparsity_out(int);
int {{ ocp.model_name }}_expl_vde_forw_n_in();
int {{ ocp.model_name }}_expl_vde_forw_n_out();

// explicit adjoint VDE
int {{ ocp.model_name }}_expl_vde_adj(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ocp.model_name }}_expl_vde_adj_work(int *, int *, int *, int *);
const int *{{ ocp.model_name }}_expl_vde_adj_sparsity_in(int);
const int *{{ ocp.model_name }}_expl_vde_adj_sparsity_out(int);
int {{ ocp.model_name }}_expl_vde_adj_n_in();
int {{ ocp.model_name }}_expl_vde_adj_n_out();

// explicit adjoint ODE jac
int {{ ocp.model_name }}_expl_ode_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ocp.model_name }}_expl_ode_hess_work(int *, int *, int *, int *);
const int *{{ ocp.model_name }}_expl_ode_hess_sparsity_in(int);
const int *{{ ocp.model_name }}_expl_ode_hess_sparsity_out(int);
int {{ ocp.model_name }}_expl_ode_hess_n_in();
int {{ ocp.model_name }}_expl_ode_hess_n_out();

{% endif %}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // {{ ocp.model_name }}_MODEL
