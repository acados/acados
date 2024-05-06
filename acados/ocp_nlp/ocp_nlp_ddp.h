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
/// \addtogroup ocp_nlp_ddp ocp_nlp_ddp
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_DDP_H_
#define ACADOS_OCP_NLP_OCP_NLP_DDP_H_

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
    double tol_zero_res; // exit tolerance if objective function is 0 for least-squares problem
    int max_iter;
    int ext_qp_res;      // compute external QP residuals (i.e. at SQP level) at each SQP iteration (for debugging)
    int qp_warm_start;   // qp_warm_start in all but the first ddp iterations
    bool warm_start_first_qp; // to set qp_warm_start in first iteration
    int rti_phase;       // only phase 0 at the moment
    int initialize_t_slacks;  // 0-false or 1-true

    // Line search
    double linesearch_eta;
    double linesearch_minimum_step_size;
    double linesearch_step_size_reduction_factor;

    // Flag for usage of adaptive levenberg marquardt strategy
    bool with_adaptive_levenberg_marquardt;
    double adaptive_levenberg_marquardt_lam;
    double adaptive_levenberg_marquardt_mu_min;
    double adaptive_levenberg_marquardt_mu0;

} ocp_nlp_ddp_opts;

//
acados_size_t ocp_nlp_ddp_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_ddp_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_ddp_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_ddp_opts_update(void *config, void *dims, void *opts);
//
void ocp_nlp_ddp_opts_set(void *config_, void *opts_, const char *field, void* value);
//
void ocp_nlp_ddp_opts_set_at_stage(void *config_, void *opts_, size_t stage, const char *field, void* value);


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

    // statistics
    double *stat;
    int stat_m;
    int stat_n;

    int status;
    int ddp_iter;

    // ddp specific memory
    double *tmp_nu_times_nx;
    struct blasfeo_dmat K_mat;

    // regularization for Levenberg-Marquardt
    double mu;
    double mu_bar;
    double step_norm;

} ocp_nlp_ddp_memory;

//
acados_size_t ocp_nlp_ddp_memory_calculate_size(void *config, void *dims, void *opts_);
//
void *ocp_nlp_ddp_memory_assign(void *config, void *dims, void *opts_, void *raw_memory);
//
void ocp_nlp_ddp_memory_reset_qp_solver(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
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

} ocp_nlp_ddp_workspace;

//
acados_size_t ocp_nlp_ddp_workspace_calculate_size(void *config, void *dims, void *opts_);



/************************************************
 * functions
 ************************************************/

//
int ocp_nlp_ddp(void *config, void *dims, void *nlp_in, void *nlp_out,
                void *args, void *mem, void *work_);
//
void ocp_nlp_ddp_config_initialize_default(void *config_);
//
int ocp_nlp_ddp_precompute(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_);
//
void ocp_nlp_ddp_eval_lagr_grad_p(void *config_, void *dims_, void *nlp_in_, void *opts_, void *mem_, void *work_,
                            const char *field, void *lagr_grad_wrt_params);
//
void ocp_nlp_ddp_get(void *config_, void *dims_, void *mem_, const char *field, void *return_value_);

int ocp_nlp_ddp_backtracking_line_search(void *config, void *dims, void *nlp_in, void *nlp_out,
                void *args, void *mem, void *work_);

double ocp_nlp_ddp_compute_qp_objective_value(ocp_nlp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out,
                ocp_nlp_workspace *nlp_work, ocp_nlp_memory *nlp_mem);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_DDP_H_
/// @}
/// @}
/// @}
