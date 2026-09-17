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


/// \addtogroup ocp_nlp
/// @{
/// \addtogroup ocp_nlp_cost ocp_nlp_cost
/// @{
/// \addtogroup ocp_nlp_cost_conl ocp_nlp_cost_conl
/// \brief This module implements convex-over-nonlinear costs of the form
/// \f$\min_{x,u,z} \psi(y(x,u,z,p) - y_{\text{ref}}, p)\f$,


#ifndef ACADOS_OCP_NLP_OCP_NLP_COST_CONL_H_
#define ACADOS_OCP_NLP_OCP_NLP_COST_CONL_H_

#ifdef __cplusplus
extern "C" {
#endif

// blasfeo
#include "blasfeo_common.h"

// acados
#include "acados/ocp_nlp/ocp_nlp_cost_common.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/types.h"



/************************************************
 * model
 ************************************************/

typedef struct
{
    // slack penalty has the form z^T * s + .5 * s^T * Z * s
    external_function_generic *conl_cost_fun;
    external_function_generic *conl_cost_fun_jac_hess;
    struct blasfeo_dvec y_ref;
    ocp_nlp_cost_common_model *common;  ///< fields shared across cost modules
    double t; // time (always zero) to match signature of external function wrt cost integration
} ocp_nlp_cost_conl_model;

//
acados_size_t ocp_nlp_cost_conl_model_calculate_size(void *config, void *dims);
//
void *ocp_nlp_cost_conl_model_assign(void *config, void *dims, void *raw_memory);
//
int ocp_nlp_cost_conl_model_set(void *config_, void *dims_, void *model_, const char *field, void *value_);



/************************************************
 * options
 ************************************************/

typedef struct
{
    bool gauss_newton_hess;  // dummy options, we always use a gauss-newton hessian
    int integrator_cost; // > 0 indicating that cost is propagated within integrator instead of cost module, only add slack contributions
    int add_hess_contribution;
} ocp_nlp_cost_conl_opts;

//
acados_size_t ocp_nlp_cost_conl_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_cost_conl_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_cost_conl_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_cost_conl_opts_update(void *config, void *dims, void *opts);
//
void ocp_nlp_cost_conl_opts_set(void *config, void *opts, const char *field, void *value);



/************************************************
 * memory
 ************************************************/
typedef struct
{
    ocp_nlp_cost_common_memory *common;  ///< fields shared across cost modules
    struct blasfeo_dmat W_chol;        // cholesky factor of hessian of outer loss function
    struct blasfeo_dvec W_chol_diag;   // cholesky factor of hessian of outer loss function if Hessian is diagonal
        // NOTE: could be in work, but needed for compatibility with NLS and cost integration
    double outer_hess_is_diag;
} ocp_nlp_cost_conl_memory;

//
acados_size_t ocp_nlp_cost_conl_memory_calculate_size(void *config, void *dims, void *opts);
//
void *ocp_nlp_cost_conl_memory_assign(void *config, void *dims, void *opts, void *raw_memory);
//
void *ocp_nlp_cost_conl_memory_get(void *memory_, const char *field);

/************************************************
 * workspace
 ************************************************/

typedef struct
{
    struct blasfeo_dmat W;             // hessian of outer loss function
    struct blasfeo_dmat Jt_ux;         // jacobian of inner residual function
    struct blasfeo_dmat Jt_ux_tilde;   // jacobian of inner residual function plus gradient contribution of algebraic variables
    struct blasfeo_dmat Jt_z;          // jacobian of inner residual function wrt algebraic variables
    struct blasfeo_dmat tmp_nv_ny;
    struct blasfeo_dmat tmp_nv_ny2;
    struct blasfeo_dmat J_y_tilde;     // workspace for integrator cost
    struct blasfeo_dvec tmp_ny;
    struct blasfeo_dvec tmp_2ns;
} ocp_nlp_cost_conl_workspace;

//
acados_size_t ocp_nlp_cost_conl_workspace_calculate_size(void *config, void *dims, void *opts);
//
size_t ocp_nlp_cost_conl_get_external_fun_workspace_requirement(void *config_, void *dims_, void *opts_, void *model_);
//
void ocp_nlp_cost_conl_set_external_fun_workspaces(void *config_, void *dims_, void *opts_, void *model_, void *workspace_);


/************************************************
 * functions
 ************************************************/

//
void ocp_nlp_cost_conl_precompute(void *config_, void *dims_, void *model_, void *opts_, void *memory_, void *work_);
//
void ocp_nlp_cost_conl_config_initialize_default(void *config, int stage);
//
void ocp_nlp_cost_conl_initialize(void *config_, void *dims, void *model_, void *opts_, void *mem_, void *work_);
//
void ocp_nlp_cost_conl_update_qp_matrices(void *config_, void *dims, void *model_, void *opts_, void *memory_, void *work_);
//
void ocp_nlp_cost_conl_compute_fun(void *config_, void *dims, void *model_, void *opts_, void *memory_, void *work_);
//
void ocp_nlp_cost_conl_compute_jac_p(void *config_, void *dims, void *model_, void *opts_, void *memory_, void *work_);
//
void ocp_nlp_cost_conl_eval_grad_p(void *config_, void *dims, void *model_, void *opts_, void *memory_, void *work_, struct blasfeo_dvec *out);
//
// Adds the contribution of one collocation node of one integrator step to the
// integrator cost. Accumulates: does NOT zero cost_fun/cost_grad/cost_hess.
void ocp_nlp_cost_conl_add_integrator_stage_cost(void *cost_capsule,
        struct blasfeo_dvec *xt, struct blasfeo_dvec *u, struct blasfeo_dvec_args *z_alg,
        struct blasfeo_dmat *S_forw_stage, double t_current, double weight,
        struct blasfeo_dmat *cost_hess);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_COST_CONL_H_
/// @}
/// @}
/// @}
