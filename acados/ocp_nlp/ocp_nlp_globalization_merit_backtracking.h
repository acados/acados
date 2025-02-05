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
/// \addtogroup ocp_nlp_globalization
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_GLOBALIZATION_MERIT_BACKTRACKING_H_
#define ACADOS_OCP_NLP_OCP_NLP_GLOBALIZATION_MERIT_BACKTRACKING_H_

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
// TODO: remove stufF!
typedef struct
{
    ocp_nlp_globalization_opts *globalization_opts;

} ocp_nlp_globalization_merit_backtracking_opts;

//
acados_size_t ocp_nlp_globalization_merit_backtracking_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_globalization_merit_backtracking_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_globalization_merit_backtracking_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_globalization_merit_backtracking_opts_set(void *config, void *opts, const char *field, void* value);


/************************************************
 * memory
 ************************************************/

// TODO: remove stufF!
typedef struct
{
    double step_norm;
    double alpha;
} ocp_nlp_globalization_merit_backtracking_memory;

//
acados_size_t ocp_nlp_globalization_merit_backtracking_memory_calculate_size(void *config, void *dims);
//
void *ocp_nlp_globalization_merit_backtracking_memory_assign(void *config, void *dims, void *raw_memory);
//

/************************************************
 * functions
 ************************************************/
//
double ocp_nlp_compute_merit_gradient(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                  ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
                                  ocp_nlp_memory *mem, ocp_nlp_workspace *work);
//
int ocp_nlp_line_search(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work,
            int sqp_iter, double *alpha_ref);
//
double ocp_nlp_evaluate_merit_fun(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
          ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work);
//
void merit_backtracking_initialize_weights(ocp_nlp_dims *dims, ocp_nlp_out *weight_merit_fun, ocp_qp_out *qp_out);
//
void merit_backtracking_update_weights(ocp_nlp_dims *dims, ocp_nlp_out *weight_merit_fun, ocp_qp_out *qp_out);
//
int ocp_nlp_globalization_merit_backtracking_find_acceptable_iterate(void *nlp_config_, void *nlp_dims_, void *nlp_in_, void *nlp_out_, void *nlp_mem_, void *solver_mem, void *nlp_work_, void *nlp_opts_, double *step_size);
//
int ocp_nlp_globalization_merit_backtracking_find_acceptable_iterate_for_ddp(void *nlp_config_, void *nlp_dims_, void *nlp_in_, void *nlp_out_, void *nlp_mem_, void *solver_mem, void *nlp_work_, void *nlp_opts_, double *step_size);
//
void ocp_nlp_globalization_merit_backtracking_print_iteration_header();
//
void ocp_nlp_globalization_merit_backtracking_print_iteration(double objective_value,
                                                            int iter_count,
                                                            void *nlp_res_,
                                                            double step_norm,
                                                            double reg_param,
                                                            int qp_status,
                                                            int qp_iter,
                                                            void* nlp_opts_,
                                                            void* mem_);
//
int ocp_nlp_globalization_merit_backtracking_needs_objective_value();
//
int ocp_nlp_globalization_merit_backtracking_needs_qp_objective_value();
//
int ocp_nlp_globalization_merit_backtracking_ddp_needs_qp_objective_value();
//
void ocp_nlp_globalization_merit_backtracking_initialize_memory(void *config_,
    void *dims_, void *nlp_mem_, void *nlp_opts_);
//
void ocp_nlp_globalization_merit_backtracking_config_initialize_default(ocp_nlp_globalization_config *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_GLOBALIZATION_MERIT_BACKTRACKING_H_
/// @}
/// @}
