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
    void (*print_iteration)(double objective_value,
                            int iter_count,
                            void* nlp_res_,
                            double step_norm,
                            double reg_param,
                            int qp_status,
                            int qp_iter,
                            void *globalization_opts,
                            void* globalization_mem);
    int (*needs_objective_value)();
    int (*needs_qp_objective_value)();
    void (*initialize_memory)(void *config_,
                            void *dims_,
                            void *nlp_mem_,
                            void *nlp_opts_);
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
