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


///
/// \defgroup ocp_nlp_cost ocp_nlp_cost
///

/// \addtogroup ocp_nlp_cost ocp_nlp_cost
/// @{
/// \addtogroup ocp_nlp_cost_common ocp_nlp_cost_common
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_COST_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_COST_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

// blasfeo
#include "blasfeo_common.h"

// acados
#include "acados/utils/external_function_generic.h"
#include "acados/utils/types.h"



/************************************************
 * dims
 ************************************************/
typedef struct
{
    int nx;  // number of states
    int nz;  // number of algebraic variables
    int nu;  // number of inputs
    int ny;  // number of outputs
    int ns;  // number of slacks
    int np;
    int np_global;
} ocp_nlp_cost_dims;


acados_size_t ocp_nlp_cost_dims_calculate_size(void *config);

void *ocp_nlp_cost_dims_assign(void *config, void *raw_memory);
//
void ocp_nlp_cost_dims_set(void *config_, void *dims_, const char *field, int* value);
//
void ocp_nlp_cost_dims_get(void *config_, void *dims_, const char *field, int* value);



/************************************************
 * common model
 ************************************************/

/// structure containing model fields shared across cost modules
typedef struct
{
    struct blasfeo_dvec Z_usr;          ///< user-provided diagonal Hessian of slacks (lower and upper)
    struct blasfeo_dvec z_usr;          ///< user-provided gradient of slacks (lower and upper)
    struct blasfeo_dvec Z_nlp;          ///< NLP-adjusted diagonal Hessian of slacks (lower and upper)
    struct blasfeo_dvec z_nlp;          ///< NLP-adjusted gradient of slacks (lower and upper)
    double scaling;                     ///< cost scaling factor
} ocp_nlp_cost_common_model;

//
acados_size_t ocp_nlp_cost_common_model_calculate_size(ocp_nlp_cost_dims* dims);
//
ocp_nlp_cost_common_model *ocp_nlp_cost_common_model_assign(ocp_nlp_cost_dims* dims, char **c_ptr);
//
int ocp_nlp_cost_common_model_set(ocp_nlp_cost_dims *dims, ocp_nlp_cost_common_model *model, const char *field, void *value_);
//
int ocp_nlp_cost_common_model_get(ocp_nlp_cost_dims *dims, ocp_nlp_cost_common_model *model, const char *field, void *value_);
//


/************************************************
 * common memory
 ************************************************/

/// structure containing memory fields shared across cost modules
typedef struct
{
    struct blasfeo_dvec grad;           ///< gradient of cost function
    struct blasfeo_dvec *ux;            ///< pointer to ux in nlp_out
    struct blasfeo_dvec *z_alg;         ///< pointer to z in sim_out
    struct blasfeo_dmat *dzdux_tran;    ///< pointer to sensitivity of z wrt ux in sim_out
    struct blasfeo_dmat *RSQrq;         ///< pointer to RSQrq in qp_in
    struct blasfeo_dvec *Z;             ///< pointer to Z in qp_in
    struct blasfeo_dvec *orphan_mask;   ///< pointer to orphan_mask in NLP memory
    struct blasfeo_dvec *seed_ux;    // pointer
    struct blasfeo_dmat *jac_lag_stat_p_global;    // pointer to jacobian of stationarity condition wrt parameters
    struct blasfeo_dvec *adj_lag_p_global;    // pointer to OCP adjoint wrt parameters
    double fun;                         ///< value of the cost function
} ocp_nlp_cost_common_memory;

//
acados_size_t ocp_nlp_cost_common_memory_calculate_size(ocp_nlp_cost_dims *dims);
//
ocp_nlp_cost_common_memory *ocp_nlp_cost_common_memory_assign(ocp_nlp_cost_dims *dims, char **c_ptr);
//
double *ocp_nlp_cost_common_memory_get_fun_ptr(ocp_nlp_cost_common_memory *memory);
//
struct blasfeo_dvec *ocp_nlp_cost_common_memory_get_grad_ptr(ocp_nlp_cost_common_memory *memory);
//
int ocp_nlp_cost_common_memory_set(ocp_nlp_cost_common_memory *memory, const char *field, void *value);


/************************************************
 * options
 ************************************************/

/// structure containing the options shared across cost modules
typedef struct
{
    int compute_hess;                   ///< if > 0, compute a Hessian approximation of cost (can only be turned off for LLS)
    int exact_hess;                     ///< if > 0, compute exact Hessian instead of Gauss-Newton approximation
    int use_numerical_hessian;          ///< > 0 indicating custom Hessian is used instead of CasADi evaluation (external cost)
    int integrator_cost;                ///< > 0 indicating that cost is propagated within integrator instead of cost module, only add slack contributions
    int with_solution_sens_wrt_params_forw;  ///< > 0 indicating that forward solution sensitivities wrt params can be computed (external cost)
    int with_solution_sens_wrt_params_adj;   ///< > 0 indicating that adjoint solution sensitivities wrt params can be computed (external cost)
    int add_hess_contribution;          ///< if > 0, add Hessian contribution to existing Hessian instead of overwriting it
} ocp_nlp_cost_common_opts;

//
acados_size_t ocp_nlp_cost_common_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_cost_common_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_cost_common_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_cost_common_opts_set(void *config, void *opts, const char *field, void *value);
//
int *ocp_nlp_cost_common_opts_get_add_hess_contribution_ptr(void *config, void *opts);



/************************************************
 * config
 ************************************************/

typedef struct
{
    acados_size_t (*dims_calculate_size)(void *config);
    void *(*dims_assign)(void *config, void *raw_memory);
    void (*dims_set)(void *config_, void *dims_, const char *field, int *value);
    void (*dims_get)(void *config_, void *dims_, const char *field, int *value);
    acados_size_t (*model_calculate_size)(void *config, void *dims);
    void *(*model_assign)(void *config, void *dims, void *raw_memory);
    int (*model_set)(void *config_, void *dims_, void *model_, const char *field, void *value_);
    int (*model_get)(void *config_, void *dims_, void *model_, const char *field, void *value_);
    acados_size_t (*opts_calculate_size)(void *config, void *dims);
    void *(*opts_assign)(void *config, void *dims, void *raw_memory);
    void (*opts_initialize_default)(void *config, void *dims, void *opts);
    void (*opts_update)(void *config, void *dims, void *opts);
    void (*opts_set)(void *config, void *opts, const char *field, void *value);
    int *(*opts_get_add_hess_contribution_ptr)(void *config, void *opts);
    acados_size_t (*memory_calculate_size)(void *config, void *dims, void *opts);
    void *(*memory_get)(void *memory_, const char *field);
    struct blasfeo_dvec *(*model_get_y_ref_ptr)(void *memory);
    double *(*model_get_scaling_ptr)(void *memory);
    double *(*get_outer_hess_is_diag_ptr)(void *memory_, void *model_);
    void (*memory_set)(void *config_, void *dims_, void *memory_, const char *field, void *value);
    void *(*memory_assign)(void *config, void *dims, void *opts, void *raw_memory);
    acados_size_t (*workspace_calculate_size)(void *config, void *dims, void *opts);
    acados_size_t (*get_external_fun_workspace_requirement)(void *config, void *dims, void *opts_, void *in);
    void (*set_external_fun_workspaces)(void *config, void *dims, void *opts_, void *in, void *work_);
    void (*initialize)(void *config_, void *dims, void *model_, void *opts_, void *mem_, void *work_);

    // computes the function value, gradient and hessian (approximation) of the cost function
    void (*update_qp_matrices)(void *config_, void *dims, void *model_, void *opts_, void *mem_, void *work_);
    // computes the cost function value (intended for globalization)
    void (*compute_fun)(void *config_, void *dims, void *model_, void *opts_, void *mem_, void *work_);
    // computes the cost jacobian wrt parameters (intended for solution sensitivities)
    void (*compute_jac_p)(void *config_, void *dims, void *model_, void *opts_, void *mem_, void *work_);
    void (*compute_adj_sol_sens_pdiff)(void *config_, void *dims, void *model_, void *opts_, void *mem_, void *work_);
    void (*eval_grad_p)(void *config_, void *dim, void* model, void *opts, void *mem, void *work, struct blasfeo_dvec *out);
    void (*compute_gradient)(void *config_, void *dims, void *model_, void *opts_, void *mem_, void *work_);
    void (*config_initialize_default)(void *config, int stage);
    void (*precompute)(void *config_, void *dims_, void *model_, void *opts_, void *memory_, void *work_);
    // stage information
    int stage;
} ocp_nlp_cost_config;

//
acados_size_t ocp_nlp_cost_config_calculate_size();
//
ocp_nlp_cost_config *ocp_nlp_cost_config_assign(void *raw_memory);


/* common functionality */
void ocp_nlp_cost_common_initialize(ocp_nlp_cost_dims *dims, ocp_nlp_cost_common_model *model, ocp_nlp_cost_common_memory *memory);
void cost_common_add_slack_contributions_to_fun_and_scale(ocp_nlp_cost_dims *dims, ocp_nlp_cost_common_model *model, ocp_nlp_cost_common_memory *memory, struct blasfeo_dvec *tmp_2ns);
void cost_common_update_gradient_with_slacks_and_scale(ocp_nlp_cost_dims *dims, ocp_nlp_cost_common_model *model, ocp_nlp_cost_common_memory *memory);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_COST_COMMON_H_
/// @}
/// @}
