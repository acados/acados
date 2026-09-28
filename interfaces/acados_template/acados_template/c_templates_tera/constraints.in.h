/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */

#ifndef {{ model.name }}_CONSTRAINTS
#define {{ model.name }}_CONSTRAINTS

#ifdef __cplusplus
extern "C" {
#endif

{% if dims.nphi > 0 %}
int {{ model.name }}_phi_constraint_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_phi_constraint_fun_work(int *, int *, int *, int *);
const int *{{ model.name }}_phi_constraint_fun_sparsity_in(int);
const int *{{ model.name }}_phi_constraint_fun_sparsity_out(int);
int {{ model.name }}_phi_constraint_fun_n_in(void);
int {{ model.name }}_phi_constraint_fun_n_out(void);

int {{ model.name }}_phi_constraint_fun_jac_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_phi_constraint_fun_jac_hess_work(int *, int *, int *, int *);
const int *{{ model.name }}_phi_constraint_fun_jac_hess_sparsity_in(int);
const int *{{ model.name }}_phi_constraint_fun_jac_hess_sparsity_out(int);
int {{ model.name }}_phi_constraint_fun_jac_hess_n_in(void);
int {{ model.name }}_phi_constraint_fun_jac_hess_n_out(void);
{% endif %}

{% if dims.nphi_e > 0 %}
int {{ model.name }}_phi_e_constraint_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_phi_e_constraint_fun_work(int *, int *, int *, int *);
const int *{{ model.name }}_phi_e_constraint_fun_sparsity_in(int);
const int *{{ model.name }}_phi_e_constraint_fun_sparsity_out(int);
int {{ model.name }}_phi_e_constraint_fun_n_in(void);
int {{ model.name }}_phi_e_constraint_fun_n_out(void);

int {{ model.name }}_phi_e_constraint_fun_jac_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_phi_e_constraint_fun_jac_hess_work(int *, int *, int *, int *);
const int *{{ model.name }}_phi_e_constraint_fun_jac_hess_sparsity_in(int);
const int *{{ model.name }}_phi_e_constraint_fun_jac_hess_sparsity_out(int);
int {{ model.name }}_phi_e_constraint_fun_jac_hess_n_in(void);
int {{ model.name }}_phi_e_constraint_fun_jac_hess_n_out(void);
{% endif %}

{% if dims.nphi_0 > 0 %}
int {{ model.name }}_phi_0_constraint_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_phi_0_constraint_fun_work(int *, int *, int *, int *);
const int *{{ model.name }}_phi_0_constraint_fun_sparsity_in(int);
const int *{{ model.name }}_phi_0_constraint_fun_sparsity_out(int);
int {{ model.name }}_phi_0_constraint_fun_n_in(void);
int {{ model.name }}_phi_0_constraint_fun_n_out(void);

int {{ model.name }}_phi_0_constraint_fun_jac_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_phi_0_constraint_fun_jac_hess_work(int *, int *, int *, int *);
const int *{{ model.name }}_phi_0_constraint_fun_jac_hess_sparsity_in(int);
const int *{{ model.name }}_phi_0_constraint_fun_jac_hess_sparsity_out(int);
int {{ model.name }}_phi_0_constraint_fun_jac_hess_n_in(void);
int {{ model.name }}_phi_0_constraint_fun_jac_hess_n_out(void);
{% endif %}


{% if dims.nh > 0 %}
int {{ model.name }}_constr_h_fun_jac_uxt_zt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_fun_jac_uxt_zt_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_fun_jac_uxt_zt_sparsity_in(int);
const int *{{ model.name }}_constr_h_fun_jac_uxt_zt_sparsity_out(int);
int {{ model.name }}_constr_h_fun_jac_uxt_zt_n_in(void);
int {{ model.name }}_constr_h_fun_jac_uxt_zt_n_out(void);

int {{ model.name }}_constr_h_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_fun_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_fun_sparsity_in(int);
const int *{{ model.name }}_constr_h_fun_sparsity_out(int);
int {{ model.name }}_constr_h_fun_n_in(void);
int {{ model.name }}_constr_h_fun_n_out(void);

{% if code_gen_options.with_solution_sens_wrt_params_forw %}
int {{ model.name }}_constr_h_jac_p_hess_xu_p(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_jac_p_hess_xu_p_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_jac_p_hess_xu_p_sparsity_in(int);
const int *{{ model.name }}_constr_h_jac_p_hess_xu_p_sparsity_out(int);
int {{ model.name }}_constr_h_jac_p_hess_xu_p_n_in(void);
int {{ model.name }}_constr_h_jac_p_hess_xu_p_n_out(void);
{% endif %}

{% if code_gen_options.with_solution_sens_wrt_params_adj %}
int {{ model.name }}_constr_h_hess_ux_pdiff_adj_pdiff(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_hess_ux_pdiff_adj_pdiff_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_hess_ux_pdiff_adj_pdiff_sparsity_in(int);
const int *{{ model.name }}_constr_h_hess_ux_pdiff_adj_pdiff_sparsity_out(int);
int {{ model.name }}_constr_h_hess_ux_pdiff_adj_pdiff_n_in(void);
int {{ model.name }}_constr_h_hess_ux_pdiff_adj_pdiff_n_out(void);
{% endif %}

{% if code_gen_options.with_value_sens_wrt_params %}
int {{ model.name }}_constr_h_adj_p(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_adj_p_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_adj_p_sparsity_in(int);
const int *{{ model.name }}_constr_h_adj_p_sparsity_out(int);
int {{ model.name }}_constr_h_adj_p_n_in(void);
int {{ model.name }}_constr_h_adj_p_n_out(void);
{% endif %}

{% if solver_options.hessian_approx == "EXACT" -%}
int {{ model.name }}_constr_h_fun_jac_uxt_zt_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_fun_jac_uxt_zt_hess_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_fun_jac_uxt_zt_hess_sparsity_in(int);
const int *{{ model.name }}_constr_h_fun_jac_uxt_zt_hess_sparsity_out(int);
int {{ model.name }}_constr_h_fun_jac_uxt_zt_hess_n_in(void);
int {{ model.name }}_constr_h_fun_jac_uxt_zt_hess_n_out(void);
{% endif %}
{% endif %}

{% if dims.nh_0 > 0 %}
int {{ model.name }}_constr_h_0_fun_jac_uxt_zt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_0_fun_jac_uxt_zt_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_0_fun_jac_uxt_zt_sparsity_in(int);
const int *{{ model.name }}_constr_h_0_fun_jac_uxt_zt_sparsity_out(int);
int {{ model.name }}_constr_h_0_fun_jac_uxt_zt_n_in(void);
int {{ model.name }}_constr_h_0_fun_jac_uxt_zt_n_out(void);

int {{ model.name }}_constr_h_0_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_0_fun_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_0_fun_sparsity_in(int);
const int *{{ model.name }}_constr_h_0_fun_sparsity_out(int);
int {{ model.name }}_constr_h_0_fun_n_in(void);
int {{ model.name }}_constr_h_0_fun_n_out(void);

{% if code_gen_options.with_solution_sens_wrt_params_forw %}
int {{ model.name }}_constr_h_0_jac_p_hess_xu_p(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_0_jac_p_hess_xu_p_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_0_jac_p_hess_xu_p_sparsity_in(int);
const int *{{ model.name }}_constr_h_0_jac_p_hess_xu_p_sparsity_out(int);
int {{ model.name }}_constr_h_0_jac_p_hess_xu_p_n_in(void);
int {{ model.name }}_constr_h_0_jac_p_hess_xu_p_n_out(void);
{% endif %}

{% if code_gen_options.with_solution_sens_wrt_params_adj %}
int {{ model.name }}_constr_h_0_hess_ux_pdiff_adj_pdiff(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_0_hess_ux_pdiff_adj_pdiff_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_0_hess_ux_pdiff_adj_pdiff_sparsity_in(int);
const int *{{ model.name }}_constr_h_0_hess_ux_pdiff_adj_pdiff_sparsity_out(int);
int {{ model.name }}_constr_h_0_hess_ux_pdiff_adj_pdiff_n_in(void);
int {{ model.name }}_constr_h_0_hess_ux_pdiff_adj_pdiff_n_out(void);
{% endif %}


{% if code_gen_options.with_value_sens_wrt_params %}
int {{ model.name }}_constr_h_0_adj_p(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_0_adj_p_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_0_adj_p_sparsity_in(int);
const int *{{ model.name }}_constr_h_0_adj_p_sparsity_out(int);
int {{ model.name }}_constr_h_0_adj_p_n_in(void);
int {{ model.name }}_constr_h_0_adj_p_n_out(void);
{% endif %}

{% if solver_options.hessian_approx == "EXACT" -%}
int {{ model.name }}_constr_h_0_fun_jac_uxt_zt_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_0_fun_jac_uxt_zt_hess_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_0_fun_jac_uxt_zt_hess_sparsity_in(int);
const int *{{ model.name }}_constr_h_0_fun_jac_uxt_zt_hess_sparsity_out(int);
int {{ model.name }}_constr_h_0_fun_jac_uxt_zt_hess_n_in(void);
int {{ model.name }}_constr_h_0_fun_jac_uxt_zt_hess_n_out(void);
{% endif %}
{% endif %}


{% if dims.nh_e > 0 %}
int {{ model.name }}_constr_h_e_fun_jac_uxt_zt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_e_fun_jac_uxt_zt_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_e_fun_jac_uxt_zt_sparsity_in(int);
const int *{{ model.name }}_constr_h_e_fun_jac_uxt_zt_sparsity_out(int);
int {{ model.name }}_constr_h_e_fun_jac_uxt_zt_n_in(void);
int {{ model.name }}_constr_h_e_fun_jac_uxt_zt_n_out(void);

int {{ model.name }}_constr_h_e_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_e_fun_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_e_fun_sparsity_in(int);
const int *{{ model.name }}_constr_h_e_fun_sparsity_out(int);
int {{ model.name }}_constr_h_e_fun_n_in(void);
int {{ model.name }}_constr_h_e_fun_n_out(void);

{% if code_gen_options.with_solution_sens_wrt_params_forw %}
int {{ model.name }}_constr_h_e_jac_p_hess_xu_p(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_e_jac_p_hess_xu_p_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_e_jac_p_hess_xu_p_sparsity_in(int);
const int *{{ model.name }}_constr_h_e_jac_p_hess_xu_p_sparsity_out(int);
int {{ model.name }}_constr_h_e_jac_p_hess_xu_p_n_in(void);
int {{ model.name }}_constr_h_e_jac_p_hess_xu_p_n_out(void);
{% endif %}

{% if code_gen_options.with_solution_sens_wrt_params_adj %}
int {{ model.name }}_constr_h_e_hess_ux_pdiff_adj_pdiff(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_e_hess_ux_pdiff_adj_pdiff_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_e_hess_ux_pdiff_adj_pdiff_sparsity_in(int);
const int *{{ model.name }}_constr_h_e_hess_ux_pdiff_adj_pdiff_sparsity_out(int);
int {{ model.name }}_constr_h_e_hess_ux_pdiff_adj_pdiff_n_in(void);
int {{ model.name }}_constr_h_e_hess_ux_pdiff_adj_pdiff_n_out(void);
{% endif %}


{% if code_gen_options.with_value_sens_wrt_params %}
int {{ model.name }}_constr_h_e_adj_p(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_e_adj_p_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_e_adj_p_sparsity_in(int);
const int *{{ model.name }}_constr_h_e_adj_p_sparsity_out(int);
int {{ model.name }}_constr_h_e_adj_p_n_in(void);
int {{ model.name }}_constr_h_e_adj_p_n_out(void);
{% endif %}

{% if solver_options.hessian_approx == "EXACT" -%}
int {{ model.name }}_constr_h_e_fun_jac_uxt_zt_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_constr_h_e_fun_jac_uxt_zt_hess_work(int *, int *, int *, int *);
const int *{{ model.name }}_constr_h_e_fun_jac_uxt_zt_hess_sparsity_in(int);
const int *{{ model.name }}_constr_h_e_fun_jac_uxt_zt_hess_sparsity_out(int);
int {{ model.name }}_constr_h_e_fun_jac_uxt_zt_hess_n_in(void);
int {{ model.name }}_constr_h_e_fun_jac_uxt_zt_hess_n_out(void);
{% endif %}
{% endif %}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // {{ model.name }}_CONSTRAINTS
