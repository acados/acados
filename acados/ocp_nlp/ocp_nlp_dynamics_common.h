/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


/// \ingroup ocp_nlp
/// @{

/// \defgroup ocp_nlp_dynamics ocp_nlp_dynamics
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif



// blasfeo
#include "blasfeo_common.h"

// acados
#include "acados/sim/sim_common.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/types.h"



/************************************************
 * config
 ************************************************/

typedef struct
{
    void (*config_initialize_default)(void *config, int stage);
    sim_config *sim_solver;
    /* dims */
    acados_size_t (*dims_calculate_size)(void *config);
    void *(*dims_assign)(void *config, void *raw_memory);
    void (*dims_set)(void *config_, void *dims_, const char *field, int *value);
    void (*dims_get)(void *config_, void *dims_, const char *field, int* value);
    /* model */
    acados_size_t (*model_calculate_size)(void *config, void *dims);
    void *(*model_assign)(void *config, void *dims, void *raw_memory);
    void (*model_set)(void *config_, void *dims_, void *model_, const char *field, void *value_);
    /* opts */
    acados_size_t (*opts_calculate_size)(void *config, void *dims);
    void *(*opts_assign)(void *config, void *dims, void *raw_memory);
    void (*opts_initialize_default)(void *config, void *dims, void *opts);
    void (*opts_set)(void *config_, void *opts_, const char *field, void *value);
    void (*opts_get)(void *config_, void *opts_, const char *field, void *value);
    void (*opts_update)(void *config, void *dims, void *opts);
    /* memory */
    acados_size_t (*memory_calculate_size)(void *config, void *dims, void *opts);
    void *(*memory_assign)(void *config, void *dims, void *opts, void *raw_memory);
    // get shooting node gap x_next(x_n, u_n) - x_{n+1}
    struct blasfeo_dvec *(*memory_get_fun_ptr)(void *memory_);
    struct blasfeo_dvec *(*memory_get_adj_ptr)(void *memory_);
    void (*memory_get)(void *config, void *dims, void *mem, const char *field, void* value);
    void (*memory_set)(void *config, void *dims, void *mem, const char *field, void* value);
    /* workspace */
    acados_size_t (*workspace_calculate_size)(void *config, void *dims, void *opts);
    void (*initialize)(void *config_, void *dims, void *model_, void *opts_, void *mem_, void *work_);
    void (*update_qp_matrices)(void *config_, void *dims, void *model_, void *opts_, void *mem_, void *work_);
    void (*compute_fun)(void *config_, void *dims, void *model_, void *opts_, void *mem_, void *work_);
    void (*compute_jac_hess_p)(void *config_, void *dims, void *model_, void *opts, void *mem, void *work_);
    void (*compute_adj_sol_sens_pdiff)(void *config_, void *dims, void *model_, void *opts, void *mem, void *work_);

    acados_size_t (*get_external_fun_workspace_requirement)(void *config, void *dims, void *opts_, void *in);
    void (*set_external_fun_workspaces)(void *config, void *dims, void *opts_, void *in, void *work_);

    void (*compute_adj_p)(void *config, void *dims, void *model, void *opts, void *memory, struct blasfeo_dvec *out);
    void (*compute_fun_and_adj)(void *config_, void *dims, void *model_, void *opts_, void *mem_, void *work_);
    int (*precompute)(void *config_, void *dims, void *model_, void *opts_, void *mem_, void *work_);
    void (*reset)(void *config_, void *dims_, void *model_, void *opts_, void *mem_, void *work_);
    int stage;
} ocp_nlp_dynamics_config;

//
acados_size_t ocp_nlp_dynamics_config_calculate_size();
//
ocp_nlp_dynamics_config *ocp_nlp_dynamics_config_assign(void *raw_memory);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_COMMON_H_
/// @}
/// @}
