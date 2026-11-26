/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef {{ model.name }}_COST
#define {{ model.name }}_COST

#ifdef __cplusplus
extern "C" {
#endif

{% if model.has_slack_cost_expr.initial %}
int {{ model.name }}_cost_slack_penalty_fun_zZ_0(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_slack_penalty_fun_zZ_0_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_slack_penalty_fun_zZ_0_sparsity_in(int);
const int *{{ model.name }}_cost_slack_penalty_fun_zZ_0_sparsity_out(int);
int {{ model.name }}_cost_slack_penalty_fun_zZ_0_n_in(void);
int {{ model.name }}_cost_slack_penalty_fun_zZ_0_n_out(void);
{% if code_gen_options.with_solution_sens_wrt_params_forw %}
int {{ model.name }}_cost_slack_penalty_hess_s_p_0(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_slack_penalty_hess_s_p_0_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_slack_penalty_hess_s_p_0_sparsity_in(int);
const int *{{ model.name }}_cost_slack_penalty_hess_s_p_0_sparsity_out(int);
int {{ model.name }}_cost_slack_penalty_hess_s_p_0_n_in(void);
int {{ model.name }}_cost_slack_penalty_hess_s_p_0_n_out(void);
{% endif %}
{% if code_gen_options.with_solution_sens_wrt_params_adj %}
int {{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_0(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_0_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_0_sparsity_in(int);
const int *{{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_0_sparsity_out(int);
int {{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_0_n_in(void);
int {{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_0_n_out(void);
{% endif %}
{% if code_gen_options.with_value_sens_wrt_params %}
int {{ model.name }}_cost_slack_penalty_fun_grad_p_0(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_slack_penalty_fun_grad_p_0_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_slack_penalty_fun_grad_p_0_sparsity_in(int);
const int *{{ model.name }}_cost_slack_penalty_fun_grad_p_0_sparsity_out(int);
int {{ model.name }}_cost_slack_penalty_fun_grad_p_0_n_in(void);
int {{ model.name }}_cost_slack_penalty_fun_grad_p_0_n_out(void);
{% endif %}
{% endif %}
{% if model.has_slack_cost_expr.path %}
int {{ model.name }}_cost_slack_penalty_fun_zZ(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_slack_penalty_fun_zZ_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_slack_penalty_fun_zZ_sparsity_in(int);
const int *{{ model.name }}_cost_slack_penalty_fun_zZ_sparsity_out(int);
int {{ model.name }}_cost_slack_penalty_fun_zZ_n_in(void);
int {{ model.name }}_cost_slack_penalty_fun_zZ_n_out(void);
{% if code_gen_options.with_solution_sens_wrt_params_forw %}
int {{ model.name }}_cost_slack_penalty_hess_s_p(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_slack_penalty_hess_s_p_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_slack_penalty_hess_s_p_sparsity_in(int);
const int *{{ model.name }}_cost_slack_penalty_hess_s_p_sparsity_out(int);
int {{ model.name }}_cost_slack_penalty_hess_s_p_n_in(void);
int {{ model.name }}_cost_slack_penalty_hess_s_p_n_out(void);
{% endif %}
{% if code_gen_options.with_solution_sens_wrt_params_adj %}
int {{ model.name }}_cost_slack_penalty_fun_zZ_adj_p(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_sparsity_in(int);
const int *{{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_sparsity_out(int);
int {{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_n_in(void);
int {{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_n_out(void);
{% endif %}
{% if code_gen_options.with_value_sens_wrt_params %}
int {{ model.name }}_cost_slack_penalty_fun_grad_p(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_slack_penalty_fun_grad_p_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_slack_penalty_fun_grad_p_sparsity_in(int);
const int *{{ model.name }}_cost_slack_penalty_fun_grad_p_sparsity_out(int);
int {{ model.name }}_cost_slack_penalty_fun_grad_p_n_in(void);
int {{ model.name }}_cost_slack_penalty_fun_grad_p_n_out(void);
{% endif %}
{% endif %}
{% if model.has_slack_cost_expr.terminal %}
int {{ model.name }}_cost_slack_penalty_fun_zZ_e(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_slack_penalty_fun_zZ_e_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_slack_penalty_fun_zZ_e_sparsity_in(int);
const int *{{ model.name }}_cost_slack_penalty_fun_zZ_e_sparsity_out(int);
int {{ model.name }}_cost_slack_penalty_fun_zZ_e_n_in(void);
int {{ model.name }}_cost_slack_penalty_fun_zZ_e_n_out(void);
{% if code_gen_options.with_solution_sens_wrt_params_forw %}
int {{ model.name }}_cost_slack_penalty_hess_s_p_e(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_slack_penalty_hess_s_p_e_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_slack_penalty_hess_s_p_e_sparsity_in(int);
const int *{{ model.name }}_cost_slack_penalty_hess_s_p_e_sparsity_out(int);
int {{ model.name }}_cost_slack_penalty_hess_s_p_e_n_in(void);
int {{ model.name }}_cost_slack_penalty_hess_s_p_e_n_out(void);
{% endif %}
{% if code_gen_options.with_solution_sens_wrt_params_adj %}
int {{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_e(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_e_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_e_sparsity_in(int);
const int *{{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_e_sparsity_out(int);
int {{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_e_n_in(void);
int {{ model.name }}_cost_slack_penalty_fun_zZ_adj_p_e_n_out(void);
{% endif %}
{% if code_gen_options.with_value_sens_wrt_params %}
int {{ model.name }}_cost_slack_penalty_fun_grad_p_e(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_slack_penalty_fun_grad_p_e_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_slack_penalty_fun_grad_p_e_sparsity_in(int);
const int *{{ model.name }}_cost_slack_penalty_fun_grad_p_e_sparsity_out(int);
int {{ model.name }}_cost_slack_penalty_fun_grad_p_e_n_in(void);
int {{ model.name }}_cost_slack_penalty_fun_grad_p_e_n_out(void);
{% endif %}
{% endif %}

// Cost at initial shooting node
{% if cost.cost_type_0 == "NONLINEAR_LS" %}
int {{ model.name }}_cost_y_0_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_y_0_fun_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_y_0_fun_sparsity_in(int);
const int *{{ model.name }}_cost_y_0_fun_sparsity_out(int);
int {{ model.name }}_cost_y_0_fun_n_in(void);
int {{ model.name }}_cost_y_0_fun_n_out(void);

int {{ model.name }}_cost_y_0_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_y_0_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_y_0_fun_jac_ut_xt_sparsity_in(int);
const int *{{ model.name }}_cost_y_0_fun_jac_ut_xt_sparsity_out(int);
int {{ model.name }}_cost_y_0_fun_jac_ut_xt_n_in(void);
int {{ model.name }}_cost_y_0_fun_jac_ut_xt_n_out(void);

{%- if solver_options.hessian_approx == "EXACT" %}
int {{ model.name }}_cost_y_0_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_y_0_hess_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_y_0_hess_sparsity_in(int);
const int *{{ model.name }}_cost_y_0_hess_sparsity_out(int);
int {{ model.name }}_cost_y_0_hess_n_in(void);
int {{ model.name }}_cost_y_0_hess_n_out(void);
{%- endif %}
{% elif cost.cost_type_0 == "CONVEX_OVER_NONLINEAR" %}

int {{ model.name }}_conl_cost_0_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_conl_cost_0_fun_work(int *, int *, int *, int *);
const int *{{ model.name }}_conl_cost_0_fun_sparsity_in(int);
const int *{{ model.name }}_conl_cost_0_fun_sparsity_out(int);
int {{ model.name }}_conl_cost_0_fun_n_in(void);
int {{ model.name }}_conl_cost_0_fun_n_out(void);

int {{ model.name }}_conl_cost_0_fun_jac_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_conl_cost_0_fun_jac_hess_work(int *, int *, int *, int *);
const int *{{ model.name }}_conl_cost_0_fun_jac_hess_sparsity_in(int);
const int *{{ model.name }}_conl_cost_0_fun_jac_hess_sparsity_out(int);
int {{ model.name }}_conl_cost_0_fun_jac_hess_n_in(void);
int {{ model.name }}_conl_cost_0_fun_jac_hess_n_out(void);

{% elif cost.cost_type_0 == "EXTERNAL" %}
    {%- if cost.cost_ext_fun_type_0 == "casadi" %}
int {{ model.name }}_cost_ext_cost_0_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_0_fun_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_0_fun_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_0_fun_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_0_fun_n_in(void);
int {{ model.name }}_cost_ext_cost_0_fun_n_out(void);

int {{ model.name }}_cost_ext_cost_0_fun_jac_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_0_fun_jac_hess_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_0_fun_jac_hess_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_0_fun_jac_hess_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_0_fun_jac_hess_n_in(void);
int {{ model.name }}_cost_ext_cost_0_fun_jac_hess_n_out(void);

int {{ model.name }}_cost_ext_cost_0_fun_jac(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_0_fun_jac_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_0_fun_jac_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_0_fun_jac_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_0_fun_jac_n_in(void);
int {{ model.name }}_cost_ext_cost_0_fun_jac_n_out(void);

{% if code_gen_options.with_solution_sens_wrt_params_forw %}
int {{ model.name }}_cost_ext_cost_0_hess_xu_p(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_0_hess_xu_p_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_0_hess_xu_p_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_0_hess_xu_p_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_0_hess_xu_p_n_in(void);
int {{ model.name }}_cost_ext_cost_0_hess_xu_p_n_out(void);
{% endif %}

{% if code_gen_options.with_solution_sens_wrt_params_adj %}
int {{ model.name }}_cost_ext_cost_0_adj_ux_pdiff(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_0_adj_ux_pdiff_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_0_adj_ux_pdiff_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_0_adj_ux_pdiff_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_0_adj_ux_pdiff_n_in(void);
int {{ model.name }}_cost_ext_cost_0_adj_ux_pdiff_n_out(void);
{% endif %}

{% if code_gen_options.with_value_sens_wrt_params %}
int {{ model.name }}_cost_ext_cost_0_grad_p(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_0_grad_p_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_0_grad_p_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_0_grad_p_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_0_grad_p_n_in(void);
int {{ model.name }}_cost_ext_cost_0_grad_p_n_out(void);
{% endif %}

    {%- else %}
int {{ cost.cost_generic_fun_0 }}(void **, void **, void *);
int {{ cost.cost_generic_fun_jac_0 }}(void **, void **, void *);
{%- if solver_options.hessian_approx == "EXACT" %}
int {{ cost.cost_generic_fun_jac_hess_0 }}(void **, void **, void *);
{%- endif %}

    {%- endif %}
{% endif %}


// Cost at path shooting node
{% if cost.cost_type == "NONLINEAR_LS" %}
int {{ model.name }}_cost_y_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_y_fun_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_y_fun_sparsity_in(int);
const int *{{ model.name }}_cost_y_fun_sparsity_out(int);
int {{ model.name }}_cost_y_fun_n_in(void);
int {{ model.name }}_cost_y_fun_n_out(void);

int {{ model.name }}_cost_y_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_y_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_y_fun_jac_ut_xt_sparsity_in(int);
const int *{{ model.name }}_cost_y_fun_jac_ut_xt_sparsity_out(int);
int {{ model.name }}_cost_y_fun_jac_ut_xt_n_in(void);
int {{ model.name }}_cost_y_fun_jac_ut_xt_n_out(void);

{%- if solver_options.hessian_approx == "EXACT" %}
int {{ model.name }}_cost_y_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_y_hess_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_y_hess_sparsity_in(int);
const int *{{ model.name }}_cost_y_hess_sparsity_out(int);
int {{ model.name }}_cost_y_hess_n_in(void);
int {{ model.name }}_cost_y_hess_n_out(void);
{%- endif %}

{% elif cost.cost_type == "CONVEX_OVER_NONLINEAR" %}
int {{ model.name }}_conl_cost_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_conl_cost_fun_work(int *, int *, int *, int *);
const int *{{ model.name }}_conl_cost_fun_sparsity_in(int);
const int *{{ model.name }}_conl_cost_fun_sparsity_out(int);
int {{ model.name }}_conl_cost_fun_n_in(void);
int {{ model.name }}_conl_cost_fun_n_out(void);

int {{ model.name }}_conl_cost_fun_jac_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_conl_cost_fun_jac_hess_work(int *, int *, int *, int *);
const int *{{ model.name }}_conl_cost_fun_jac_hess_sparsity_in(int);
const int *{{ model.name }}_conl_cost_fun_jac_hess_sparsity_out(int);
int {{ model.name }}_conl_cost_fun_jac_hess_n_in(void);
int {{ model.name }}_conl_cost_fun_jac_hess_n_out(void);
{% elif cost.cost_type == "EXTERNAL" %}
    {%- if cost.cost_ext_fun_type == "casadi" %}
int {{ model.name }}_cost_ext_cost_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_fun_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_fun_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_fun_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_fun_n_in(void);
int {{ model.name }}_cost_ext_cost_fun_n_out(void);

int {{ model.name }}_cost_ext_cost_fun_jac_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_fun_jac_hess_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_fun_jac_hess_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_fun_jac_hess_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_fun_jac_hess_n_in(void);
int {{ model.name }}_cost_ext_cost_fun_jac_hess_n_out(void);

int {{ model.name }}_cost_ext_cost_fun_jac(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_fun_jac_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_fun_jac_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_fun_jac_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_fun_jac_n_in(void);
int {{ model.name }}_cost_ext_cost_fun_jac_n_out(void);

{% if code_gen_options.with_solution_sens_wrt_params_forw %}
int {{ model.name }}_cost_ext_cost_hess_xu_p(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_hess_xu_p_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_hess_xu_p_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_hess_xu_p_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_hess_xu_p_n_in(void);
int {{ model.name }}_cost_ext_cost_hess_xu_p_n_out(void);
{% endif %}

{% if code_gen_options.with_solution_sens_wrt_params_adj %}
int {{ model.name }}_cost_ext_cost_adj_ux_pdiff(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_adj_ux_pdiff_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_adj_ux_pdiff_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_adj_ux_pdiff_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_adj_ux_pdiff_n_in(void);
int {{ model.name }}_cost_ext_cost_adj_ux_pdiff_n_out(void);
{% endif %}

{% if code_gen_options.with_value_sens_wrt_params %}
int {{ model.name }}_cost_ext_cost_grad_p(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_grad_p_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_grad_p_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_grad_p_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_grad_p_n_in(void);
int {{ model.name }}_cost_ext_cost_grad_p_n_out(void);
{% endif %}

    {%- else %}
int {{ cost.cost_generic_fun }}(void **, void **, void *);
int {{ cost.cost_generic_fun_jac }}(void **, void **, void *);
{%- if solver_options.hessian_approx == "EXACT" %}
int {{ cost.cost_generic_fun_jac_hess }}(void **, void **, void *);
{%- endif %}
    {%- endif %}
{% endif %}

// Cost at terminal shooting node
{% if cost.cost_type_e == "NONLINEAR_LS" %}
int {{ model.name }}_cost_y_e_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_y_e_fun_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_y_e_fun_sparsity_in(int);
const int *{{ model.name }}_cost_y_e_fun_sparsity_out(int);
int {{ model.name }}_cost_y_e_fun_n_in(void);
int {{ model.name }}_cost_y_e_fun_n_out(void);

int {{ model.name }}_cost_y_e_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_y_e_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_y_e_fun_jac_ut_xt_sparsity_in(int);
const int *{{ model.name }}_cost_y_e_fun_jac_ut_xt_sparsity_out(int);
int {{ model.name }}_cost_y_e_fun_jac_ut_xt_n_in(void);
int {{ model.name }}_cost_y_e_fun_jac_ut_xt_n_out(void);

{%- if solver_options.hessian_approx == "EXACT" %}
int {{ model.name }}_cost_y_e_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_y_e_hess_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_y_e_hess_sparsity_in(int);
const int *{{ model.name }}_cost_y_e_hess_sparsity_out(int);
int {{ model.name }}_cost_y_e_hess_n_in(void);
int {{ model.name }}_cost_y_e_hess_n_out(void);
{%- endif %}

{% elif cost.cost_type_e == "CONVEX_OVER_NONLINEAR" %}
int {{ model.name }}_conl_cost_e_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_conl_cost_e_fun_work(int *, int *, int *, int *);
const int *{{ model.name }}_conl_cost_e_fun_sparsity_in(int);
const int *{{ model.name }}_conl_cost_e_fun_sparsity_out(int);
int {{ model.name }}_conl_cost_e_fun_n_in(void);
int {{ model.name }}_conl_cost_e_fun_n_out(void);

int {{ model.name }}_conl_cost_e_fun_jac_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_conl_cost_e_fun_jac_hess_work(int *, int *, int *, int *);
const int *{{ model.name }}_conl_cost_e_fun_jac_hess_sparsity_in(int);
const int *{{ model.name }}_conl_cost_e_fun_jac_hess_sparsity_out(int);
int {{ model.name }}_conl_cost_e_fun_jac_hess_n_in(void);
int {{ model.name }}_conl_cost_e_fun_jac_hess_n_out(void);
{% elif cost.cost_type_e == "EXTERNAL" %}
    {%- if cost.cost_ext_fun_type_e == "casadi" %}
int {{ model.name }}_cost_ext_cost_e_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_e_fun_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_e_fun_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_e_fun_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_e_fun_n_in(void);
int {{ model.name }}_cost_ext_cost_e_fun_n_out(void);

int {{ model.name }}_cost_ext_cost_e_fun_jac_hess(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_e_fun_jac_hess_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_e_fun_jac_hess_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_e_fun_jac_hess_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_e_fun_jac_hess_n_in(void);
int {{ model.name }}_cost_ext_cost_e_fun_jac_hess_n_out(void);

int {{ model.name }}_cost_ext_cost_e_fun_jac(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_e_fun_jac_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_e_fun_jac_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_e_fun_jac_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_e_fun_jac_n_in(void);
int {{ model.name }}_cost_ext_cost_e_fun_jac_n_out(void);

{% if code_gen_options.with_solution_sens_wrt_params_forw %}
int {{ model.name }}_cost_ext_cost_e_hess_xu_p(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_e_hess_xu_p_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_e_hess_xu_p_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_e_hess_xu_p_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_e_hess_xu_p_n_in(void);
int {{ model.name }}_cost_ext_cost_e_hess_xu_p_n_out(void);
{% endif %}

{% if code_gen_options.with_solution_sens_wrt_params_adj %}
int {{ model.name }}_cost_ext_cost_e_adj_ux_pdiff(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_e_adj_ux_pdiff_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_e_adj_ux_pdiff_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_e_adj_ux_pdiff_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_e_adj_ux_pdiff_n_in(void);
int {{ model.name }}_cost_ext_cost_e_adj_ux_pdiff_n_out(void);
{% endif %}

{% if code_gen_options.with_value_sens_wrt_params %}
int {{ model.name }}_cost_ext_cost_e_grad_p(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ model.name }}_cost_ext_cost_e_grad_p_work(int *, int *, int *, int *);
const int *{{ model.name }}_cost_ext_cost_e_grad_p_sparsity_in(int);
const int *{{ model.name }}_cost_ext_cost_e_grad_p_sparsity_out(int);
int {{ model.name }}_cost_ext_cost_e_grad_p_n_in(void);
int {{ model.name }}_cost_ext_cost_e_grad_p_n_out(void);
{% endif %}

    {%- else %}
int {{ cost.cost_generic_fun_e }}(void **, void **, void *);
int {{ cost.cost_generic_fun_jac_e }}(void **, void **, void *);
{%- if solver_options.hessian_approx == "EXACT" %}
int {{ cost.cost_generic_fun_jac_hess_e }}(void **, void **, void *);
{%- endif %}
    {%- endif %}
{% endif %}


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // {{ model.name }}_COST
