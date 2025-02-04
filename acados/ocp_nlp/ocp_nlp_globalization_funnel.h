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

#ifndef ACADOS_OCP_NLP_OCP_NLP_GLOBALIZATION_FUNNEL_H_
#define ACADOS_OCP_NLP_OCP_NLP_GLOBALIZATION_FUNNEL_H_

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
    ocp_nlp_globalization_opts *globalization_opts;

    // Funnel globalization related options
    double initialization_increase_factor; // for multiplication with initial infeasibility
    double initialization_upper_bound; // for initialization of initial funnel width
    double sufficient_decrease_factor; // multiplication factor for funnel suff. decrease factor
    double kappa; // parameter for reduction of funnel
    double fraction_switching_condition; // parameter in switching condition
    double initial_penalty_parameter; // initial penalty parameter for penalty phase
    double penalty_eta; // fraction in penalty update
    double penalty_contraction; // penalty contraction factor
    bool type_switching_condition; // which type of switching condition do we use?
} ocp_nlp_globalization_funnel_opts;

//
acados_size_t ocp_nlp_globalization_funnel_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_globalization_funnel_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_globalization_funnel_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_globalization_funnel_opts_set(void *config, void *opts, const char *field, void* value);


/************************************************
 * memory
 ************************************************/

typedef struct
{
    double funnel_width;
    char funnel_iter_type;
    bool funnel_penalty_mode;
    double l1_infeasibility;
    double penalty_parameter;
    double alpha;

} ocp_nlp_globalization_funnel_memory;
//
acados_size_t ocp_nlp_globalization_funnel_memory_calculate_size(void *config, void *dims);
//
void *ocp_nlp_globalization_funnel_memory_assign(void *config, void *dims, void *raw_memory);
//
/************************************************
 * functions
 ************************************************/

void debug_output(ocp_nlp_opts *opts, char* message, int print_level);
//
void debug_output_double(ocp_nlp_opts *opts, char* message, double value, int print_level);
//
void initialize_funnel_width(ocp_nlp_globalization_funnel_memory *mem, ocp_nlp_globalization_funnel_opts *opts, double initial_infeasibility);
//
void update_funnel_penalty_parameter(ocp_nlp_globalization_funnel_memory *mem,
                                            ocp_nlp_globalization_funnel_opts *opts,
                                            double pred_f, double pred_h);
//
void decrease_funnel(ocp_nlp_globalization_funnel_memory *mem, ocp_nlp_globalization_funnel_opts *opts, double trial_infeasibility, double current_infeasibility);
//
bool is_iterate_inside_of_funnel(ocp_nlp_globalization_funnel_memory *mem, ocp_nlp_globalization_funnel_opts *opts, double infeasibility);
//
bool is_funnel_sufficient_decrease_satisfied(ocp_nlp_globalization_funnel_memory *mem, ocp_nlp_globalization_funnel_opts *opts, double infeasibility);
//
bool is_switching_condition_satisfied(ocp_nlp_globalization_funnel_opts *opts, double pred_optimality, double step_size, double pred_infeasibility);
//
bool is_f_type_armijo_condition_satisfied(ocp_nlp_globalization_opts *globalization_opts,
                                        double negative_ared,
                                        double pred,
                                        double alpha);
//
bool is_trial_iterate_acceptable_to_funnel(ocp_nlp_globalization_funnel_memory *mem,
                                            ocp_nlp_opts *nlp_opts,
                                            double pred, double ared, double alpha,
                                            double current_infeasibility,
                                            double trial_infeasibility,
                                            double current_objective,
                                            double trial_objective,
                                            double current_merit,
                                            double trial_merit,
                                            double pred_merit);
//
int backtracking_line_search(ocp_nlp_config *config,
                            ocp_nlp_dims *dims,
                            ocp_nlp_in *nlp_in,
                            ocp_nlp_out *nlp_out,
                            ocp_nlp_memory *nlp_mem,
                            void *solver_mem,
                            ocp_nlp_workspace *nlp_work,
                            ocp_nlp_opts *nlp_opts,
                            double *step_size);
//
int ocp_nlp_globalization_funnel_find_acceptable_iterate(void *nlp_config_, void *nlp_dims_, void *nlp_in_, void *nlp_out_, void *nlp_mem_, void *solver_mem, void *nlp_work_, void *nlp_opts_, double *step_size);
//
void ocp_nlp_globalization_funnel_print_iteration_header();
//
void ocp_nlp_globalization_funnel_print_iteration(double objective_value,
                                                int iter_count,
                                                void* nlp_res_,
                                                double step_norm,
                                                double reg_param,
                                                int qp_status,
                                                int qp_iter,
                                                void* nlp_opts_,
                                                void* mem_);
//
int ocp_nlp_globalization_funnel_needs_objective_value();
//
int ocp_nlp_globalization_funnel_needs_qp_objective_value();
//
void ocp_nlp_globalization_funnel_initialize_memory(void *config_, void *dims_, void *nlp_mem_, void *nlp_opts_);
//
void ocp_nlp_globalization_funnel_config_initialize_default(ocp_nlp_globalization_config *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_GLOBALIZATION_FUNNEL_H_
/// @}
/// @}
