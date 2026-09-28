/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */



/// \defgroup ocp_nlp ocp_nlp
/// @{

/// \defgroup ocp_nlp_globalization ocp_nlp_globalization
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_GLOBALIZATION_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_GLOBALIZATION_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

// This would cause cyclic include is not possible due to cycle
// #include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/utils/types.h"

/************************************************
 * config
 ************************************************/

typedef struct
{
    /* opts */
    acados_size_t (*opts_calculate_size)(void *config, void *dims);
    void *(*opts_assign)(void *config, void *dims, void *raw_memory);
    void (*opts_initialize_default)(void *config, void *dims, void *opts);
    void (*opts_set)(void *config, void *opts, const char *field, void* value);
    /* memory */
    acados_size_t (*memory_calculate_size)(void *config, void *dims);
    void *(*memory_assign)(void *config, void *dims, void *raw_memory);
    /* functions */
    int (*find_acceptable_iterate)(void *nlp_config, void *nlp_dims, void *nlp_in, void *nlp_out, void *nlp_mem, void *solver_mem, void *nlp_work, void *nlp_opts, double *step_size);
    void (*print_iteration_header)();
    void (*print_iteration)(double objective_value, void *globalization_opts, void* globalization_mem);
    int (*needs_objective_value)();
    int (*needs_qp_objective_value)();
    void (*initialize_memory)(void *config_, void *dims_, void *nlp_mem_, void *nlp_opts_);
} ocp_nlp_globalization_config;

//
acados_size_t ocp_nlp_globalization_config_calculate_size();
//
ocp_nlp_globalization_config *ocp_nlp_globalization_config_assign(void *raw_memory);


/************************************************
 * options
 ************************************************/
typedef struct ocp_nlp_globalization_opts
{
    int use_SOC;
    int line_search_use_sufficient_descent;
    int full_step_dual;
    double alpha_min;
    double alpha_reduction;
    double eps_sufficient_descent;
} ocp_nlp_globalization_opts;

//
acados_size_t ocp_nlp_globalization_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_globalization_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_globalization_opts_initialize_default(void *config, void *dims, void *opts);
//
// void ocp_nlp_globalization_opts_update(void *config, void *dims, void *opts);
//
void ocp_nlp_globalization_opts_set(void *config_, void *opts_, const char *field, void* value);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_GLOBALIZATION_COMMON_H_
/// @}
/// @}
