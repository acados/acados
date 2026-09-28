/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


/// \ingroup ocp_nlp
/// @{

/// \defgroup ocp_nlp_constraints ocp_nlp_constraints
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/types.h"



/************************************************
 * config
 ************************************************/

typedef struct
{
    acados_size_t (*dims_calculate_size)(void *config);
    void *(*dims_assign)(void *config, void *raw_memory);
    acados_size_t (*model_calculate_size)(void *config, void *dims);
    void *(*model_assign)(void *config, void *dims, void *raw_memory);
    int (*model_set)(void *config_, void *dims_, void *model_, const char *field, void *value);
    void (*model_get)(void *config_, void *dims_, void *model_, const char *field, void *value);
    void (*model_set_dmask_ptr)(struct blasfeo_dvec *dmask, void *model_);
    acados_size_t (*opts_calculate_size)(void *config, void *dims);
    void *(*opts_assign)(void *config, void *dims, void *raw_memory);
    void (*opts_initialize_default)(void *config, void *dims, void *opts);
    void (*opts_update)(void *config, void *dims, void *opts);
    void (*opts_set)(void *config, void *opts, char *field, void *value);
    acados_size_t (*memory_calculate_size)(void *config, void *dims, void *opts);
    struct blasfeo_dvec *(*memory_get_fun_ptr)(void *memory);
    struct blasfeo_dvec *(*memory_get_adj_ptr)(void *memory);
    void (*memory_set)(void *config_, void *dims_, void *memory_, const char *field, void *value);

    void *(*memory_assign)(void *config, void *dims, void *opts, void *raw_memory);
    acados_size_t (*workspace_calculate_size)(void *config, void *dims, void *opts);
    acados_size_t (*get_external_fun_workspace_requirement)(void *config, void *dims, void *opts_, void *in);
    void (*set_external_fun_workspaces)(void *config, void *dims, void *opts_, void *in, void *work_);
    void (*initialize)(void *config, void *dims, void *model, void *opts, void *mem, void *work);
    void (*update_slack_masks_wrt_orphans)(void *config, void *dims, void *model, void *opts, void *mem, void *work);
    //
    void (*precompute)(void *config, void *dims, void *model, void *opts, void *mem, void *work);
    void (*update_qp_matrices)(void *config, void *dims, void *model, void *opts, void *mem, void *work);
    void (*update_qp_vectors)(void *config, void *dims, void *model, void *opts, void *mem, void *work);
    void (*compute_fun)(void *config, void *dims, void *model, void *opts, void *mem, void *work);
    void (*compute_jac_hess_p)(void *config, void *dims, void *model, void *opts, void *mem, void *work);
    void (*compute_adj_p)(void *config, void *dims, void *model, void *opts, void *memory, void *work, struct blasfeo_dvec *out);
    void (*compute_adj_sol_sens_pdiff)(void *config_, void *dims, void *model_, void *opts, void *mem, void *work_);
    void (*config_initialize_default)(void *config, int stage);
    // dimension setters
    void (*dims_set)(void *config_, void *dims_, const char *field, const int *value);
    void (*dims_get)(void *config_, void *dims_, const char *field, int* value);
    // stage information
    int stage;
} ocp_nlp_constraints_config;

//
acados_size_t ocp_nlp_constraints_config_calculate_size();
//
ocp_nlp_constraints_config *ocp_nlp_constraints_config_assign(void *raw_memory);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_COMMON_H_
/// @}
/// @}
