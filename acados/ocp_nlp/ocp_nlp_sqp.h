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
/// \addtogroup ocp_nlp_solver
/// @{
/// \addtogroup ocp_nlp_sqp ocp_nlp_sqp
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_SQP_H_
#define ACADOS_OCP_NLP_OCP_NLP_SQP_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/utils/types.h"



/************************************************
 * options
 ************************************************/

typedef struct
{
    ocp_nlp_opts *nlp_opts;
    double tol_stat;     // exit tolerance on stationarity condition
    double tol_eq;       // exit tolerance on equality constraints
    double tol_ineq;     // exit tolerance on inequality constraints
    double tol_comp;     // exit tolerance on complementarity condition
    double tol_unbounded; // exit threshold when objective function seems to be unbounded
    double tol_min_step_norm; // exit tolerance for small step
    int max_iter;
    int ext_qp_res;      // compute external QP residuals (i.e. at SQP level) at each SQP iteration (for debugging)
    int log_primal_step_norm; // compute and log the max norm of the primal steps
    int qp_warm_start;   // qp_warm_start in all but the first sqp iterations
    bool warm_start_first_qp; // to set qp_warm_start in first iteration
    int rti_phase;       // only phase 0 at the moment
    bool eval_residual_at_max_iter; // if convergence should be checked after last iterations or only throw max_iter reached

    // Funnel globalization related options
    double funnel_initialization_increase_factor; // for multiplication with initial infeasibility
    double funnel_initialization_upper_bound; // for initialization of initial funnel width
    double funnel_sufficient_decrease_factor; // multiplication factor for funnel suff. decrease factor
    double funnel_kappa; // parameter for reduction of funnel
    double funnel_fraction_switching_condition; // parameter in switching condition
    double funnel_initial_penalty_parameter; // initial penalty parameter for penalty phase
    double funnel_penalty_eta; // fraction in penalty update
    double funnel_penalty_contraction; // penalty contraction factor
    bool funnel_type_switching_condition;

} ocp_nlp_sqp_opts;

//
acados_size_t ocp_nlp_sqp_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_sqp_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_sqp_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_sqp_opts_update(void *config, void *dims, void *opts);
//
void ocp_nlp_sqp_opts_set(void *config_, void *opts_, const char *field, void* value);
//
void ocp_nlp_sqp_opts_set_at_stage(void *config_, void *opts_, size_t stage, const char *field, void* value);



/************************************************
 * memory
 ************************************************/

typedef struct
{
    // nlp memory
    ocp_nlp_memory *nlp_mem;

    double time_qp_sol;
    double time_qp_solver_call;
    double time_qp_xcond;
    double time_lin;
    double time_reg;
    double time_tot;
    double time_glob;
    double time_sim;
    double time_sim_la;
    double time_sim_ad;
    double time_solution_sensitivities;
    double alpha;

    double *primal_step_norm;

    // statistics
    double *stat;
    int stat_m;
    int stat_n;

    int status;
    int sqp_iter;

    double step_norm;

    double funnel_width;
    char funnel_iter_type;
    bool funnel_penalty_mode;
    double l1_infeasibility;
    double funnel_penalty_parameter;

} ocp_nlp_sqp_memory;

//
acados_size_t ocp_nlp_sqp_memory_calculate_size(void *config, void *dims, void *opts_);
//
void *ocp_nlp_sqp_memory_assign(void *config, void *dims, void *opts_, void *raw_memory);
//
void ocp_nlp_sqp_memory_reset_qp_solver(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
    void *opts_, void *mem_, void *work_);


/************************************************
 * workspace
 ************************************************/

typedef struct
{
    ocp_nlp_workspace *nlp_work;

    // qp residuals
    ocp_qp_res *qp_res;
    ocp_qp_res_ws *qp_res_ws;

} ocp_nlp_sqp_workspace;

//
acados_size_t ocp_nlp_sqp_workspace_calculate_size(void *config, void *dims, void *opts_);



/************************************************
 * functions
 ************************************************/

//
int ocp_nlp_sqp(void *config, void *dims, void *nlp_in, void *nlp_out,
                void *args, void *mem, void *work_);
//
void ocp_nlp_sqp_config_initialize_default(void *config_);
//
int ocp_nlp_sqp_precompute(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_);
//
void ocp_nlp_sqp_eval_lagr_grad_p(void *config_, void *dims_, void *nlp_in_, void *opts_, void *mem_, void *work_,
                            const char *field, void *grad_p);
//
void ocp_nlp_sqp_get(void *config_, void *dims_, void *mem_, const char *field, void *return_value_);
//
double ocp_nlp_sqp_compute_qp_objective_value(ocp_nlp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out,
                ocp_nlp_workspace *nlp_work, ocp_nlp_memory *nlp_mem, ocp_nlp_sqp_opts *opts);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_SQP_H_
/// @}
/// @}
/// @}
