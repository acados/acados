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


#include "acados/ocp_nlp/ocp_nlp_globalization_common.h"
#include "acados/ocp_nlp/ocp_nlp_globalization_funnel.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"

// external
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#if defined(ACADOS_WITH_OPENMP)
#include <omp.h>
#endif

// blasfeo
#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"
// acados
#include "acados/utils/mem.h"


/************************************************
 * options
 ************************************************/

acados_size_t ocp_nlp_globalization_funnel_opts_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_globalization_funnel_opts);

    size += ocp_nlp_globalization_opts_calculate_size(config, dims);

    return size;
}

void ocp_nlp_globalization_funnel_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_globalization_funnel_opts *opts = opts_;
    ocp_nlp_globalization_opts *globalization_opts = opts->globalization_opts;
    ocp_nlp_globalization_config *config = config_;

    ocp_nlp_globalization_opts_initialize_default(config, dims, globalization_opts);

    // funnel method opts
    opts->initialization_increase_factor = 15.0;
    opts->initialization_upper_bound = 1.0;
    opts->sufficient_decrease_factor = 0.9;
    opts->kappa = 0.9;
    opts->fraction_switching_condition = 1e-3;
    opts->initial_penalty_parameter = 1.0;
    opts->penalty_contraction = 0.5;
    opts->penalty_eta = 1e-6;
    opts->type_switching_condition = false; // use ipopt/gould type of switching

    return;
}


void ocp_nlp_globalization_funnel_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_globalization_funnel_opts *opts = opts_;
    ocp_nlp_globalization_config *config = config_;

    if (!strcmp(field, "funnel_init_increase_factor"))
    {
        double* funnel_init_increase_factor = (double *) value;
        if (*funnel_init_increase_factor <= 1.0)
        {
            printf("\nerror: ocp_nlp_globalization_funnel_opts_set: invalid value for funnel_init_increase_factor field, need double > 1, got %f.", *funnel_init_increase_factor);
            exit(1);
        }
        opts->initialization_increase_factor = *funnel_init_increase_factor;
    }
    else if (!strcmp(field, "funnel_init_upper_bound"))
    {
        double* funnel_init_upper_bound = (double *) value;
        if (*funnel_init_upper_bound <= 0.0)
        {
            printf("\nerror: ocp_nlp_globalization_funnel_opts_set: invalid value for funnel_init_upper_bound field, need double > 0, got %f.", *funnel_init_upper_bound);
            exit(1);
        }
        opts->initialization_upper_bound = *funnel_init_upper_bound;
    }
    else if (!strcmp(field, "funnel_sufficient_decrease_factor"))
    {
        double* funnel_sufficient_decrease_factor = (double *) value;
        if (*funnel_sufficient_decrease_factor <= 0.0 || *funnel_sufficient_decrease_factor >= 1.0)
        {
            printf("\nerror: ocp_nlp_globalization_funnel_opts_set: invalid value for funnel_sufficient_decrease_factor field, need double in (0,1), got %f.", *funnel_sufficient_decrease_factor);
            exit(1);
        }
        opts->sufficient_decrease_factor = *funnel_sufficient_decrease_factor;
    }
    else if (!strcmp(field, "funnel_kappa"))
    {
        double* funnel_kappa = (double *) value;
        if (*funnel_kappa <= 0.0 || *funnel_kappa >= 1.0)
        {
            printf("\nerror: ocp_nlp_globalization_funnel_opts_set: invalid value for funnel_kappa field, need double in (0,1), got %f.", *funnel_kappa);
            exit(1);
        }
        opts->kappa = *funnel_kappa;
    }
    else if (!strcmp(field, "funnel_fraction_switching_condition"))
    {
        double* funnel_fraction_switching_condition = (double *) value;
        if (*funnel_fraction_switching_condition <= 0.0 || *funnel_fraction_switching_condition >= 1.0)
        {
            printf("\nerror: ocp_nlp_globalization_funnel_opts_set: invalid value for funnel_fraction_switching_condition field, need double in (0,1), got %f.", *funnel_fraction_switching_condition);
            exit(1);
        }
        opts->fraction_switching_condition = *funnel_fraction_switching_condition;
    }
    else if (!strcmp(field, "funnel_initial_penalty_parameter"))
    {
        double* funnel_initial_penalty_parameter = (double *) value;
        if (*funnel_initial_penalty_parameter < 0.0 || *funnel_initial_penalty_parameter > 1.0)
        {
            printf("\nerror: ocp_nlp_globalization_funnel_opts_set: invalid value for funnel_initial_penalty_parameter field, need double in [0,1], got %f.", *funnel_initial_penalty_parameter);
            exit(1);
        }
        opts->initial_penalty_parameter = *funnel_initial_penalty_parameter;
    }
    else
    {
        ocp_nlp_globalization_opts_set(config, opts->globalization_opts, field, value);
    }
    return;
}

void *ocp_nlp_globalization_funnel_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_globalization_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_globalization_funnel_opts *opts = (ocp_nlp_globalization_funnel_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_globalization_funnel_opts);

    opts->globalization_opts = ocp_nlp_globalization_opts_assign(config, dims, c_ptr);
    c_ptr += ocp_nlp_globalization_opts_calculate_size(config, dims);

    assert((char *) raw_memory + ocp_nlp_globalization_funnel_opts_calculate_size(config_, dims_) >=
           c_ptr);

    return opts;
}

/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_globalization_funnel_memory_calculate_size(void *config_, void *dims_)
{
    acados_size_t size = 0;

    size += sizeof(ocp_nlp_globalization_funnel_memory);

    return size;
}

void *ocp_nlp_globalization_funnel_memory_assign(void *config_, void *dims_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // initial align
    align_char_to(8, &c_ptr);

    ocp_nlp_globalization_funnel_memory *mem = (ocp_nlp_globalization_funnel_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_globalization_funnel_memory);

    align_char_to(8, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_globalization_funnel_memory_calculate_size(config_, dims_) >= c_ptr);

    return mem;
}

/************************************************
 * funnel functions
 ************************************************/
static double compute_gradient_directional_derivative(ocp_nlp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    // Compute the QP objective function value
    double dir_der = 0.0;
    int i, nux, ns;
    int N = dims->N;
    // Sum over stages 0 to N
    for (i = 0; i <= N; i++)
    {
        nux = dims->nx[i] + dims->nu[i];
        ns = dims->ns[i];
        // Calculate g.T d
        dir_der += blasfeo_ddot(nux, &qp_out->ux[i], 0, &qp_in->rqz[i], 0);

        // Calculate gradient of slacks
        dir_der += blasfeo_ddot(2 * ns, &qp_out->ux[i], nux, &qp_in->rqz[i], nux);
    }
    return dir_der;
}

void debug_output(ocp_nlp_opts *opts, char* message, int print_level)
{
    if (opts->print_level > print_level)
    {
        printf("%s", message); //debugging output
    }
}
void debug_output_double(ocp_nlp_opts *opts, char* message, double value, int print_level)
{
    if (opts->print_level > print_level)
    {
        printf("%s: %f\n", message, value); //debugging output
    }
}

void initialize_funnel_width(ocp_nlp_globalization_funnel_memory *mem, ocp_nlp_globalization_funnel_opts *opts, double initial_infeasibility)
{
    mem->funnel_width = fmax(opts->initialization_upper_bound,
                            opts->initialization_increase_factor*initial_infeasibility);
}

void initialize_funnel_penalty_parameter(ocp_nlp_globalization_funnel_memory *mem, ocp_nlp_globalization_funnel_opts *opts)
{
    mem->penalty_parameter = opts->initial_penalty_parameter;
}

void update_funnel_penalty_parameter(ocp_nlp_globalization_funnel_memory *mem,
                                            ocp_nlp_globalization_funnel_opts *opts,
                                            double pred_f, double pred_h)
{
    if (mem->penalty_parameter * pred_f + pred_h < opts->penalty_eta * pred_h)
    {
        mem->penalty_parameter = fmin(opts->penalty_contraction * mem->penalty_parameter,
                                             ((1-opts->penalty_eta) * pred_h) / (-pred_f));
    }
    // else: do not decrease penalty parameter
}

void decrease_funnel(ocp_nlp_globalization_funnel_memory *mem, ocp_nlp_globalization_funnel_opts *opts, double trial_infeasibility, double current_infeasibility)
{
    mem->funnel_width = (1-opts->kappa) * trial_infeasibility + opts->kappa * mem->funnel_width;
}

bool is_iterate_inside_of_funnel(ocp_nlp_globalization_funnel_memory *mem, ocp_nlp_globalization_funnel_opts *opts, double infeasibility)
{
    if (infeasibility <= mem->funnel_width)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool is_funnel_sufficient_decrease_satisfied(ocp_nlp_globalization_funnel_memory *mem, ocp_nlp_globalization_funnel_opts *opts, double infeasibility)
{
    if (infeasibility <= opts->sufficient_decrease_factor* mem->funnel_width)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool is_switching_condition_satisfied(ocp_nlp_globalization_funnel_opts *opts, double pred_optimality, double step_size, double pred_infeasibility)
{
    // if (step_size * pred_optimality >= opts->fraction_switching_condition * pred_infeasibility * pred_infeasibility)
    if (step_size * pred_optimality >= opts->fraction_switching_condition * pred_infeasibility)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool is_f_type_armijo_condition_satisfied(ocp_nlp_globalization_opts *globalization_opts,
                                                    double negative_ared,
                                                    double pred,
                                                    double alpha)
{
    if (negative_ared <= fmin(globalization_opts->eps_sufficient_descent * alpha * fmax(pred, 0) + 1e-18, 0))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool is_trial_iterate_acceptable_to_funnel(ocp_nlp_globalization_funnel_memory *mem,
                                           ocp_nlp_opts *nlp_opts,
                                                  double pred, double ared, double alpha,
                                                  double current_infeasibility,
                                                  double trial_infeasibility,
                                                  double current_objective,
                                                  double trial_objective,
                                                  double current_merit,
                                                  double trial_merit,
                                                  double pred_merit)
{
    ocp_nlp_globalization_funnel_opts *opts = nlp_opts->globalization;
    ocp_nlp_globalization_opts *globalization_opts = opts->globalization_opts;
    bool accept_step = false;
    debug_output_double(nlp_opts, "current objective", current_objective, 2);
    debug_output_double(nlp_opts, "current infeasibility", current_infeasibility, 2);
    debug_output_double(nlp_opts, "trial objective", trial_objective, 2);
    debug_output_double(nlp_opts, "trial infeasibility", trial_infeasibility, 2);
    debug_output_double(nlp_opts, "pred", pred, 2);

    if(is_iterate_inside_of_funnel(mem, opts, trial_infeasibility))
    {
        debug_output(nlp_opts, "Trial iterate is INSIDE of funnel\n", 1);
        if (!mem->funnel_penalty_mode)
        {
            debug_output(nlp_opts, "Penalty Mode not active!\n", 1);
            if (is_switching_condition_satisfied(opts, pred, alpha, current_infeasibility))
            {
                debug_output(nlp_opts, "Switching condition IS satisfied!\n", 1);
                if (is_f_type_armijo_condition_satisfied(globalization_opts, -ared, pred, alpha))
                {
                    debug_output(nlp_opts, "f-type step: Armijo condition satisfied\n", 1);
                    accept_step = true;
                    mem->funnel_iter_type = 'f';
                }
                else
                {
                    debug_output(nlp_opts, "f-type step: Armijo condition NOT satisfied\n", 1);
                }

            }
            else if (is_funnel_sufficient_decrease_satisfied(mem, opts, trial_infeasibility))
            {
                debug_output(nlp_opts, "Switching condition is NOT satisfied!\n", 1);
                debug_output(nlp_opts, "h-type step: funnel suff. decrease satisfied!\n", 1);
                accept_step = true;
                mem->funnel_iter_type = 'h';
                decrease_funnel(mem, opts, trial_infeasibility, current_infeasibility);
            }
            else
            {
                debug_output(nlp_opts, "Switching condition is NOT satisfied!\n", 1);
                debug_output(nlp_opts, "Entered penalty check!\n", 1);
                //TODO move to function and test more
                if (trial_merit <= current_merit + globalization_opts->eps_sufficient_descent * alpha * pred_merit)
                {
                    debug_output(nlp_opts, "Penalty Function accepted\n", 1);
                    accept_step = true;
                    mem->funnel_iter_type = 'b';
                    mem->funnel_penalty_mode = true;
                }
            }
        }
        else
        {
            debug_output(nlp_opts, "Penalty mode active\n", 1);
            if (trial_merit <= current_merit + globalization_opts->eps_sufficient_descent * alpha * pred_merit)
            {
                debug_output(nlp_opts, "p-type step: accepted iterate\n", 1);
                accept_step = true;
                mem->funnel_iter_type = 'p';

                if (is_funnel_sufficient_decrease_satisfied(mem, opts, trial_infeasibility))
                {
                    decrease_funnel(mem, opts, trial_infeasibility, current_infeasibility);
                    mem->funnel_penalty_mode = false;
                }
            }
        }
    }
    else
    {
        debug_output(nlp_opts, "Trial iterate is NOT INSIDE of funnel\n", 1);
    }
    return accept_step;
}

int backtracking_line_search(ocp_nlp_config *config,
                            ocp_nlp_dims *dims,
                            ocp_nlp_in *nlp_in,
                            ocp_nlp_out *nlp_out,
                            ocp_nlp_memory *nlp_mem,
                            void *solver_mem,
                            ocp_nlp_workspace *nlp_work,
                            ocp_nlp_opts *nlp_opts,
                            double *step_size)
{
    ocp_nlp_globalization_funnel_opts *opts = nlp_opts->globalization;
    ocp_nlp_globalization_opts *globalization_opts = opts->globalization_opts;
    ocp_nlp_globalization_funnel_memory *mem = nlp_mem->globalization;

    int N = dims->N;
    double pred;
    if (opts->type_switching_condition)
    {
        pred = -nlp_mem->qp_cost_value;
    }
    else
    {
        pred = -compute_gradient_directional_derivative(dims, nlp_mem->qp_in, nlp_mem->qp_out);
    }
    double pred_merit = 0.0; // Calculate this here
    double alpha = 1.0;
    double trial_cost;
    double trial_infeasibility = 0.0;
    double ared;
    bool accept_step;
    double current_infeasibility = mem->l1_infeasibility;
    double current_cost = nlp_mem->cost_value;
    double current_merit = mem->penalty_parameter*current_cost + current_infeasibility;

    // do the penalty parameter update here .... might be changed later
    update_funnel_penalty_parameter(mem, opts, pred, mem->l1_infeasibility);

    int i;

    while (true)
    {
        // Calculate trial iterate: trial_iterate = current_iterate + alpha * direction
        config->step_update(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem,
                                     nlp_work, nlp_work->tmp_nlp_out, solver_mem, alpha, globalization_opts->full_step_dual);

        ///////////////////////////////////////////////////////////////////////
        // Evaluate cost function at trial iterate
        // set evaluation point to tmp_nlp_out
        ocp_nlp_set_primal_variable_pointers_in_submodules(config, dims, nlp_in, nlp_work->tmp_nlp_out, nlp_mem);
        // compute trial dynamics value
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
        for (i=0; i<N; i++)
        {
            // dynamics: Note has to be first, because cost_integration might be used.
            config->dynamics[i]->compute_fun(config->dynamics[i], dims->dynamics[i], nlp_in->dynamics[i],
                                            nlp_opts->dynamics[i], nlp_mem->dynamics[i], nlp_work->dynamics[i]);
        }
        // compute trial objective function value
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
        for (i=0; i<=N; i++)
        {
            // cost
            config->cost[i]->compute_fun(config->cost[i], dims->cost[i], nlp_in->cost[i], nlp_opts->cost[i],
                                        nlp_mem->cost[i], nlp_work->cost[i]);
        }
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
        for (i=0; i<=N; i++)
        {
            // constr
            config->constraints[i]->compute_fun(config->constraints[i], dims->constraints[i],
                                                nlp_in->constraints[i], nlp_opts->constraints[i],
                                                nlp_mem->constraints[i], nlp_work->constraints[i]);
        }
        // reset evaluation point to SQP iterate
        ocp_nlp_set_primal_variable_pointers_in_submodules(config, dims, nlp_in, nlp_out, nlp_mem);

        double *tmp_fun;
        // Calculate the trial objective and constraint violation
        trial_cost = 0.0;
        for(i=0; i<=N; i++)
        {
            tmp_fun = config->cost[i]->memory_get_fun_ptr(nlp_mem->cost[i]);
            trial_cost += *tmp_fun;
        }
        trial_infeasibility = ocp_nlp_get_l1_infeasibility(config, dims, nlp_mem);

        ///////////////////////////////////////////////////////////////////////
        // Evaluate merit function at trial point
        double trial_merit = mem->penalty_parameter*trial_cost + trial_infeasibility;
        pred_merit = mem->penalty_parameter * pred + current_infeasibility;
        ared = nlp_mem->cost_value - trial_cost;

        // Funnel globalization
        accept_step = is_trial_iterate_acceptable_to_funnel(mem, nlp_opts,
                                                            pred, ared,
                                                            alpha, current_infeasibility,
                                                            trial_infeasibility, current_cost,
                                                            trial_cost, current_merit, trial_merit,
                                                            pred_merit);

        if (accept_step)
        {
            mem->alpha = alpha;
            *step_size = alpha;
            nlp_mem->cost_value = trial_cost;
            mem->l1_infeasibility = trial_infeasibility;
            return ACADOS_SUCCESS;
        }

        if (alpha < globalization_opts->alpha_min)
        {
            printf("Funnel Linesearch: Step size gets too small. alpha = %e < alpha_min = %e Should enter penalty phase. \n", alpha, globalization_opts->alpha_min);
            return ACADOS_MINSTEP;
        }

        alpha *= globalization_opts->alpha_reduction;
    }
}

int ocp_nlp_globalization_funnel_find_acceptable_iterate(void *nlp_config_, void *nlp_dims_, void *nlp_in_, void *nlp_out_, void *nlp_mem_, void *solver_mem, void *nlp_work_, void *nlp_opts_, double *step_size)
{
    ocp_nlp_config *nlp_config = nlp_config_;
    ocp_nlp_dims *nlp_dims = nlp_dims_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = nlp_mem_;
    ocp_nlp_workspace *nlp_work = nlp_work_;
    ocp_nlp_opts *nlp_opts = nlp_opts_;

    int linesearch_success;
    linesearch_success = backtracking_line_search(nlp_config, nlp_dims, nlp_in, nlp_out, nlp_mem, solver_mem, nlp_work, nlp_opts, step_size);
    // Copy new iterate to nlp_out
    if (linesearch_success == ACADOS_SUCCESS)
    {
        // in case line search fails, we do not want to copy trial iterates!
        copy_ocp_nlp_out(nlp_dims, nlp_work->tmp_nlp_out, nlp_out);
    }
    return linesearch_success;
}


void ocp_nlp_globalization_funnel_print_iteration_header()
{
    printf("%6s | %11s | %10s | %10s | %10s | %10s | %10s | %10s | %10s | %12s | %10s | %10s | %10s | %10s\n",
    "iter.",
    "objective",
    "res_eq",
    "res_ineq",
    "res_stat",
    "res_comp",
    "alpha",
    "step_norm",
    "LM_reg.",
    "funnel width",
    "penalty",
    "qp_status",
    "qp_iter",
    "iter. type");
}


void ocp_nlp_globalization_funnel_print_iteration(double objective_value,
                                                int iter_count,
                                                void* nlp_res_,
                                                double step_norm,
                                                double reg_param,
                                                int qp_status,
                                                int qp_iter,
                                                void *nlp_opts_,
                                                void *mem_)
{
    ocp_nlp_res *nlp_res = nlp_res_;
    ocp_nlp_globalization_funnel_memory* mem = (ocp_nlp_globalization_funnel_memory*) mem_;
    if ((iter_count % 10 == 0))
    {
        ocp_nlp_globalization_funnel_print_iteration_header();
    }
    printf("%6i | %11.4e | %10.4e | %10.4e | %10.4e | %10.4e | %10.4e | %10.4e | %10.4e | %12.4e | %10.4e | %10i | %10i | %10c\n",
        iter_count,
        objective_value,
        nlp_res->inf_norm_res_eq,
        nlp_res->inf_norm_res_ineq,
        nlp_res->inf_norm_res_stat,
        nlp_res->inf_norm_res_comp,
        mem->alpha,
        step_norm,
        reg_param,
        mem->funnel_width,
        mem->penalty_parameter,
        qp_status,
        qp_iter,
        mem->funnel_iter_type);
}

int ocp_nlp_globalization_funnel_needs_objective_value()
{
    return 1;
}

int ocp_nlp_globalization_funnel_needs_qp_objective_value()
{
    return 1;
}

// TODO(David): maybe rename to initialize
void ocp_nlp_globalization_funnel_initialize_memory(void *config_, void *dims_, void *nlp_mem_, void *nlp_opts_)
{
    printf("Note: The funnel globalization is still under development.\n");
    printf("If you encouter problems or bugs, please report to the acados developers!\n");

    ocp_nlp_memory *nlp_mem = (ocp_nlp_memory *) nlp_mem_;
    ocp_nlp_opts *nlp_opts = (ocp_nlp_opts *) nlp_opts_;
    ocp_nlp_globalization_funnel_opts *opts = nlp_opts->globalization;
    ocp_nlp_globalization_funnel_memory *mem = nlp_mem->globalization;
    mem->l1_infeasibility = ocp_nlp_get_l1_infeasibility(config_, dims_, nlp_mem_);
    initialize_funnel_width(mem, opts, mem->l1_infeasibility);
    mem->funnel_iter_type = '-';
    initialize_funnel_penalty_parameter(mem, opts);
    mem->alpha = 1.0;
}

void ocp_nlp_globalization_funnel_config_initialize_default(ocp_nlp_globalization_config *config)
{
    // opts
    config->opts_calculate_size = &ocp_nlp_globalization_funnel_opts_calculate_size;
    config->opts_assign = &ocp_nlp_globalization_funnel_opts_assign;
    config->opts_initialize_default = &ocp_nlp_globalization_funnel_opts_initialize_default;
    config->opts_set = &ocp_nlp_globalization_funnel_opts_set;
    // memory
    config->memory_calculate_size = &ocp_nlp_globalization_funnel_memory_calculate_size;
    config->memory_assign = &ocp_nlp_globalization_funnel_memory_assign;
    // functions
    config->find_acceptable_iterate = &ocp_nlp_globalization_funnel_find_acceptable_iterate;
    config->print_iteration_header = &ocp_nlp_globalization_funnel_print_iteration_header;
    config->print_iteration = &ocp_nlp_globalization_funnel_print_iteration;
    config->needs_objective_value = &ocp_nlp_globalization_funnel_needs_objective_value;
    config->needs_qp_objective_value = &ocp_nlp_globalization_funnel_needs_qp_objective_value;
    config->initialize_memory = &ocp_nlp_globalization_funnel_initialize_memory;
}
