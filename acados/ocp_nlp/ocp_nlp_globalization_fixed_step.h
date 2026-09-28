/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


/// \addtogroup ocp_nlp
/// @{
/// \addtogroup ocp_nlp_globalization
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_GLOBALIZATION_FIXED_STEP_H_
#define ACADOS_OCP_NLP_OCP_NLP_GLOBALIZATION_FIXED_STEP_H_

#ifdef __cplusplus
extern "C" {
#endif

// blasfeo
#include "blasfeo_common.h"

// acados
#include "acados/ocp_nlp/ocp_nlp_globalization_common.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/utils/types.h"

/************************************************
 * options
 ************************************************/
typedef struct
{
    double step_length;
    ocp_nlp_globalization_opts *globalization_opts;

} ocp_nlp_globalization_fixed_step_opts;
//
acados_size_t ocp_nlp_globalization_fixed_step_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_globalization_fixed_step_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_globalization_fixed_step_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_globalization_fixed_step_opts_set(void *config, void *opts, const char *field, void* value);

/************************************************
 * memory
 ************************************************/
typedef struct
{
    double alpha;  // dummy value, not used, just to have non-empty struct
} ocp_nlp_globalization_fixed_step_memory;

//
acados_size_t ocp_nlp_globalization_fixed_step_memory_calculate_size(void *config, void *dims);
//
void *ocp_nlp_globalization_fixed_step_memory_assign(void *config, void *dims, void *raw_memory);
//

/************************************************
 * functions
 ************************************************/
//
int ocp_nlp_globalization_fixed_step_find_acceptable_iterate(void *nlp_config_, void *nlp_dims_, void *nlp_in_, void *nlp_out_, void *solver_mem, void *nlp_mem_, void *nlp_work_, void *nlp_opts_, double *step_size);
//
void ocp_nlp_globalization_fixed_step_print_iteration_header();
//
void ocp_nlp_globalization_fixed_step_print_iteration(double objective_value,
                                                    void* nlp_opts_,
                                                    void* mem_);
//
int ocp_nlp_globalization_fixed_step_needs_objective_value();
//
int ocp_nlp_globalization_fixed_step_needs_qp_objective_value();
//
void ocp_nlp_globalization_fixed_step_initialize_memory(void *config_, void *dims_, void *nlp_mem_, void *nlp_opts_);
//
void ocp_nlp_globalization_fixed_step_config_initialize_default(ocp_nlp_globalization_config *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_GLOBALIZATION_FIXED_STEP_H_
/// @}
/// @}
