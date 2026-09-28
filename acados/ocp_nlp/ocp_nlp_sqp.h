/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
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
    int log_primal_step_norm; // compute and log the max norm of the primal steps
    double timeout_max_time; // maximum time the solve may require before timeout is triggered. No timeout if 0.
    ocp_nlp_timeout_heuristic_t timeout_heuristic; // type of heuristic used to predict solve time of next QP
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

    double alpha;
    double *primal_step_norm;

    // statistics
    double *stat;
    int stat_m;
    int stat_n;

    double step_norm;
    double timeout_estimated_per_iteration_time;

} ocp_nlp_sqp_memory;

//
acados_size_t ocp_nlp_sqp_memory_calculate_size(void *config, void *dims, void *opts_, void *in_);
//
void *ocp_nlp_sqp_memory_assign(void *config, void *dims, void *opts_, void *in_, void *raw_memory);
//
void ocp_nlp_sqp_memory_reset_qp_solver(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
    void *opts_, void *mem_, void *work_);


/************************************************
 * workspace
 ************************************************/

typedef struct
{
    ocp_nlp_workspace *nlp_work;
} ocp_nlp_sqp_workspace;

//
acados_size_t ocp_nlp_sqp_workspace_calculate_size(void *config, void *dims, void *opts_, void *in_);



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

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_SQP_H_
/// @}
/// @}
/// @}
