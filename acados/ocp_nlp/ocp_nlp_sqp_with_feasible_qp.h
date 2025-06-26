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
/// \addtogroup ocp_nlp_sqp_wfqp ocp_nlp_sqp_wfqp
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_SQP_WITH_FEASIBLE_QP_H_
#define ACADOS_OCP_NLP_OCP_NLP_SQP_WITH_FEASIBLE_QP_H_

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
    bool log_pi_norm_inf; // compute and log the max norm of the pi multipliers
    bool log_lam_norm_inf; // compute and log the max norm of the lam multipliers
    bool use_constraint_hessian_in_feas_qp; // Either use exact Hessian or identity matrix in feasibility QP
    bool use_QP_l1_inf_from_slacks; // True: sums up the slack variable values from qp_out; False: compute manually; Should give the same result.
    int search_direction_mode; // determines how the QPs should be solved
    int watchdog_zero_slacks_max; // number of consecutive BYRD_OMOJOKUN iterations with zero slacks before switching back to NOMINAL_QP
    bool allow_direction_mode_switch_to_nominal; // if true, mode can switch from Byrd-Omojokun to nominal mode
    double feasibility_qp_hessian_scalar; // multiplication factor of feasibility QP Hessian
} ocp_nlp_sqp_wfqp_opts;


//
acados_size_t ocp_nlp_sqp_wfqp_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_sqp_wfqp_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_sqp_wfqp_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_sqp_wfqp_opts_update(void *config, void *dims, void *opts);
//
void ocp_nlp_sqp_wfqp_opts_set(void *config_, void *opts_, const char *field, void* value);
//
void ocp_nlp_sqp_wfqp_opts_set_at_stage(void *config_, void *opts_, size_t stage, const char *field, void* value);


/************************************************
 * memory
 ************************************************/

typedef struct
{
    // nlp memory
    ocp_nlp_memory *nlp_mem;

    double alpha;

    int *nns;  // number of non-slacked constraints in NLP
    int **idxns;  // indices of non-slacked constraints in NLP

    // statistics
    double *stat;
        // res_stat, res_eq, res_ineq, res_comp, qp_stat_N1, qp_iter_N1, qp_stat_F, qp_iter_F, qp_stat_N2, qp_iter_N2, alpha
    int stat_m;
    int stat_n;

    double step_norm;
    double norm_inf_pi;
    double norm_inf_lam;

    struct blasfeo_dvec *Z_cost_module;  // Z values from cost module
    struct blasfeo_dmat *RSQ_cost;
    struct blasfeo_dmat *RSQ_constr;

    double pred_l1_inf_QP;
    double l1_infeasibility;
    int search_direction_mode; // either NOMINAL_QP or BYRD_OMOJOKUN
    char* search_direction_type; // for output logging
    int watchdog_zero_slacks_counter; // counts number of consecutive BYRD_OMOJOKUN iter with slack sum == 0
    int absolute_nns; // sum of all nns[i]

    int qps_solved_in_iter;

    // QP solver with always feasible QPs
    ocp_qp_xcond_solver relaxed_qp_solver;
    ocp_qp_xcond_solver_memory *relaxed_qp_solver_mem;
    ocp_qp_xcond_solver_workspace *relaxed_qp_solver_work;
    // qp in & out
    ocp_qp_in *relaxed_qp_in;
    ocp_qp_out *relaxed_qp_out;
    void *relaxed_qpscaling_mem;
    // only pointers
    ocp_qp_in *relaxed_scaled_qp_in;
    ocp_qp_out *relaxed_scaled_qp_out;

} ocp_nlp_sqp_wfqp_memory;

//
acados_size_t ocp_nlp_sqp_wfqp_memory_calculate_size(void *config, void *dims, void *opts_, void *in_);
//
void *ocp_nlp_sqp_wfqp_memory_assign(void *config, void *dims, void *opts_, void *in_, void *raw_memory);
//
void ocp_nlp_sqp_wfqp_memory_reset_qp_solver(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
    void *opts_, void *mem_, void *work_);
//
void set_relaxed_qp_in_matrix_pointers(ocp_nlp_sqp_wfqp_memory *mem);

/************************************************
 * workspace
 ************************************************/

typedef struct
{
    ocp_nlp_workspace *nlp_work;
} ocp_nlp_sqp_wfqp_workspace;

//
acados_size_t ocp_nlp_sqp_wfqp_workspace_calculate_size(void *config, void *dims, void *opts_, void *in_);



/************************************************
 * functions
 ************************************************/
//
int ocp_nlp_sqp_wfqp(void *config, void *dims, void *nlp_in, void *nlp_out,
                void *args, void *mem, void *work_);
//
void ocp_nlp_sqp_wfqp_config_initialize_default(void *config_);
//
void ocp_nlp_sqp_wfqp_config_initialize_default_feasible_qp(void *config_);
//
int ocp_nlp_sqp_wfqp_precompute(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_);
//
void ocp_nlp_sqp_wfqp_eval_lagr_grad_p(void *config_, void *dims_, void *nlp_in_, void *opts_, void *mem_, void *work_,
                            const char *field, void *grad_p);
//
void ocp_nlp_sqp_wfqp_get(void *config_, void *dims_, void *mem_, const char *field, void *return_value_);
//
double ocp_nlp_sqp_wfqp_compute_qp_objective_value(ocp_nlp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out,
                ocp_nlp_workspace *nlp_work);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_SQP_WITH_FEASIBLE_QP_H_
/// @}
/// @}
/// @}