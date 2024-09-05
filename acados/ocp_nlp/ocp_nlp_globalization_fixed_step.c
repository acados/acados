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
#include "acados/ocp_nlp/ocp_nlp_globalization_fixed_step.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"

// TODO: copy boilerblate..
// fix imports

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
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/utils/mem.h"




int ocp_nlp_line_search(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work,
            int sqp_iter, double *alpha_reference)
{
    ocp_nlp_globalization_opts *globalization_opts = opts->globalization;
    int i, j;

    int N = dims->N;
    int *nv = dims->nv;

    double merit_fun1;
    ocp_qp_out *qp_out = mem->qp_out;

    if (globalization_opts->globalization == FIXED_STEP)
    {
        *alpha_reference = opts->step_length;
        return ACADOS_SUCCESS;
    }
    else if (globalization_opts->globalization != MERIT_BACKTRACKING)
    {
        printf("ocp_nlp_line_search: should only be called with globalization FIXED_STEP or MERIT_BACKTRACKING");
        exit(1);
    }

    /* MERIT_BACKTRACKING line search */
    // Following Leineweber1999, Section "3.5.1 Line Search Globalization"
    // TODO: check out more advanced step search Leineweber1995

    // copy out (current iterate) to work->tmp_nlp_out
    for (i = 0; i <= N; i++)
        blasfeo_dveccp(nv[i], out->ux+i, 0, work->tmp_nlp_out->ux+i, 0);
    // NOTE: copying duals not needed, as they dont enter the merit function

    // TODO: think about z here!
    // linear update of algebraic variables using state and input sensitivity
    //    if (i < N)
    //    {
    //        blasfeo_dgemv_t(nu[i]+nx[i], nz[i], alpha, mem->dzduxt+i, 0, 0, mem->qp_out->ux+i, 0, 1.0, mem->z_alg+i, 0, out->z+i, 0);
    //    }

    /* modify/initialize merit function weights (Leineweber1999 M5.1, p.89) */
    if (sqp_iter==0)
    {
        merit_backtracking_initialize_weights(dims, work->weight_merit_fun, qp_out);
    }
    else
    {
        merit_backtracking_update_weights(dims, work->weight_merit_fun, qp_out);
    }

    // TODO: why does Leineweber do full step in first SQP iter?
    // if (sqp_iter == 0)
    // {
    // }

    double merit_fun0 = ocp_nlp_evaluate_merit_fun(config, dims, in, out, opts, mem, work);

    double reduction_factor = globalization_opts->alpha_reduction;
    double max_next_merit_fun_val = merit_fun0;
    double eps_sufficient_descent = globalization_opts->eps_sufficient_descent;
    double dmerit_dy = 0.0;
    double alpha = 1.0;

    /* actual Line Search*/
    if (globalization_opts->line_search_use_sufficient_descent)
    {
        // check Armijo-type sufficient descent condition Leinweber1999 (2.35);
        dmerit_dy = ocp_nlp_compute_merit_gradient(config, dims, in, out, opts, mem, work);
        if (dmerit_dy > 0.0)
        {
            if (dmerit_dy > 1e-6 && opts->print_level > 0)
            {
                printf("\nacados line search: found dmerit_dy = %e > 0. Setting it to 0.0 instead\n", dmerit_dy);
            }
            dmerit_dy = 0.0;
        }
    }

    // From Leineweber1999: eq (3.64) -> only relevant for adaptive integrators looking at Remark 3.2.
    // "It is noteworthy that our practical implementation takes into account the potential nonsmoothness introduced by the fact that certain components of the penalty function - namely the continuity condition residuals - are evaluated only within integration tolerance."
    // double sum_pi = 0.0;
    // for (i = 0; i < N; i++)
    // {
    //     for (j = 0; j < dims->nx[i+1]; j++)
    //         sum_pi += BLASFEO_DVECEL(work->weight_merit_fun->pi+i, j);
    // }
    // double relaxed_val = 2.0 * 1e-6 * sum_pi;
    // if (abs(merit_fun0 - merit_fun1) < relaxed_val)
    // {
    //     printf("\nexiting because of relaxed_val.");
    //     break;
    // }

    for (j=0; alpha*reduction_factor > globalization_opts->alpha_min; j++)
    {
        // tmp_nlp_out = out + alpha * qp_out
        for (i = 0; i <= N; i++)
            blasfeo_daxpy(nv[i], alpha, qp_out->ux+i, 0, out->ux+i, 0, work->tmp_nlp_out->ux+i, 0);

        merit_fun1 = ocp_nlp_evaluate_merit_fun(config, dims, in, out, opts, mem, work);
        if (opts->print_level > 1)
        {
            printf("backtracking %d alpha = %f, merit_fun1 = %e, merit_fun0 %e\n", j, alpha, merit_fun1, merit_fun0);
        }

        // if (merit_fun1 < merit_fun0 && merit_fun1 > max_next_merit_fun_val)
        // {
        //     printf("\nalpha %f would be accepted without sufficient descent condition", alpha);
        // }

        max_next_merit_fun_val = merit_fun0 + eps_sufficient_descent * dmerit_dy * alpha;
        if ((merit_fun1 < max_next_merit_fun_val) && !isnan(merit_fun1) && !isinf(merit_fun1))
        {
            *alpha_reference = alpha;
            return ACADOS_SUCCESS;
        }
        else
        {
            alpha *= reduction_factor;
        }
    }

    *alpha_reference = alpha;
    if (isnan(merit_fun1) || isinf(merit_fun1))
    {
        return ACADOS_NAN_DETECTED;
    }
    else
    {
        return ACADOS_MINSTEP;
    }
}


void ocp_nlp_globalization_merit_backtracking_print_iteration_header()
{
    printf("# it\tstat\t\teq\t\tineq\t\tcomp\t\tqp_stat\tqp_iter\talpha\n");
}


// TODO: unified signature:
// -> move everything around.
// 1. residual_iter
// 2. int iter count
// 3. alpha etc. move to glob_memory. (void *)
void ocp_nlp_globalization_merit_backtracking_print_iteration(ocp_nlp_opts* opts,
                    double obj,
                    int iter_count,
                    double infeas_eq,
                    double infeas_ineq,
                    double stationarity,
                    double complementarity,
                    double alpha,
                    double step_norm,
                    double reg_param,
                    double funnel_width,
                    double penalty_parameter,
                    int qp_status,
                    int qp_iter,
                    char iter_type)
{
    ocp_nlp_globalization_opts *globalization_opts = opts->globalization;
    if ((iter_count % 10 == 0)){
        print_iteration_header(opts);
    }
    printf("%i\t%e\t%e\t%e\t%e\t%d\t%d\t%e\n",
        iter_count,
        stationarity,
        infeas_eq,
        infeas_ineq,
        complementarity,
        qp_status,
        qp_iter,
        alpha);
}

int ocp_nlp_globalization_step_needs_objective_value()
{
    return 0;
}

void ocp_nlp_globalization_fixed_step_config_initialize_default(ocp_nlp_globalization_config *config)
{
    // opts
    config->opts_calculate_size = &ocp_nlp_globalization_fixed_step_opts_calculate_size;
    config->opts_assign = &ocp_nlp_globalization_fixed_step_opts_assign;
    config->opts_initialize_default = &ocp_nlp_globalization_fixed_step_opts_initialize_default;
    config->opts_set = &ocp_nlp_globalization_fixed_step_opts_set;
    // functions
    config->find_acceptable_iterate = &ocp_nlp_globalization_fixed_step_find_acceptable_iterate;
    config->print_iteration_header = &ocp_nlp_globalization_fixed_step_print_iteration_header;
    config->print_iteration = &ocp_nlp_globalization_fixed_step_print_iteration;
    config->needs_objective_value = &ocp_nlp_globalization_fixed_step_needs_objective_value;
}

