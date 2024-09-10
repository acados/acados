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


#include "acados/ocp_nlp/ocp_nlp_globalization_merit_backtracking.h"

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
// acados
#include "acados/ocp_nlp/ocp_nlp_globalization_common.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
// #include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/utils/mem.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"


/************************************************
 * options
 ************************************************/

acados_size_t ocp_nlp_globalization_merit_backtracking_opts_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_globalization_merit_backtracking_opts);

    size += ocp_nlp_globalization_opts_calculate_size(config, dims);

    return size;
}

void ocp_nlp_globalization_merit_backtracking_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_globalization_opts_initialize_default(config_, dims_, opts_);
    return;
}


void ocp_nlp_globalization_merit_backtracking_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_globalization_merit_backtracking_opts *opts = opts_;
    ocp_nlp_globalization_config *config = config_;

    config->opts_set(config, opts->globalization_opts, field, value);

    return;
}

void *ocp_nlp_globalization_merit_backtracking_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_globalization_merit_backtracking_opts *opts = (ocp_nlp_globalization_merit_backtracking_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_globalization_merit_backtracking_opts);

    assert((char *) raw_memory + ocp_nlp_globalization_merit_backtracking_opts_calculate_size(config_, dims_) >=
           c_ptr);

    return opts;
}

/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_globalization_merit_backtracking_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    acados_size_t size = 0;

    size += sizeof(ocp_nlp_globalization_merit_backtracking_memory);

    return size;
}

void *ocp_nlp_globalization_merit_backtracking_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_globalization_merit_backtracking_opts *opts = opts_;

    char *c_ptr = (char *) raw_memory;

    // initial align
    align_char_to(8, &c_ptr);

    ocp_nlp_globalization_merit_backtracking_memory *mem = (ocp_nlp_globalization_merit_backtracking_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_globalization_merit_backtracking_memory);

    align_char_to(8, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_globalization_merit_backtracking_memory_calculate_size(config, dims, opts) >= c_ptr);

    return mem;
}

/************************************************
 * functions
 ************************************************/

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
                                        ocp_nlp_globalization_merit_backtracking_memory* mem)
                    // double obj,
                    // int iter_count,
                    // double infeas_eq,
                    // double infeas_ineq,
                    // double stationarity,
                    // double complementarity,
                    // double alpha,
                    // double step_norm,
                    // double reg_param,
                    // double funnel_width,
                    // double penalty_parameter,
                    // int qp_status,
                    // int qp_iter,
                    // char iter_type)
{
    // if ((iter_count % 10 == 0)){
    //     ocp_nlp_globalization_merit_backtracking_print_iteration_header();
    // }
    // printf("%i\t%e\t%e\t%e\t%e\t%d\t%d\t%e\n",
    //     iter_count,
    //     stationarity,
    //     infeas_eq,
    //     infeas_ineq,
    //     complementarity,
    //     qp_status,
    //     qp_iter,
    //     alpha);
}

static bool ocp_nlp_soc_line_search(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in,
            ocp_nlp_out *nlp_out, ocp_nlp_sqp_opts *opts, ocp_nlp_sqp_memory *mem, ocp_nlp_sqp_workspace *work, int sqp_iter)
{
    int ii;
    int N = dims->N;

    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;

    ocp_nlp_workspace *nlp_work = work->nlp_work;

    ocp_qp_in *qp_in = nlp_mem->qp_in;
    ocp_qp_out *qp_out = nlp_mem->qp_out;
    qp_info *qp_info_;
    // NOTE: following Waechter2006:
    // Do SOC
    // 1. if "the first trial step size alpha_k,0 has been rejected and
    // 2. if the infeasibility would have increased when accepting the previous step
    // NOTE: the "and" is interpreted as an "or" in the current implementation

    // preliminary line search
    int line_search_status = ocp_nlp_line_search_merit_check_full_step(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, sqp_iter);

    // return bool do_line_search;
    if (line_search_status == ACADOS_NAN_DETECTED)
    {
        // do line search but no SOC.
        return true;
    }
    else if (line_search_status == ACADOS_SUCCESS)
    {
        mem->alpha = 1.0;
        return false;
    }
    // else perform SOC (below)

    // Second Order Correction (SOC): following Nocedal2006: p.557, eq. (18.51) -- (18.56)
    // Paragraph: APPROACH III: S l1 QP (SEQUENTIAL l1 QUADRATIC PROGRAMMING),
    // Section 18.8 TRUST-REGION SQP METHODS
    //   - just no trust region radius here.
    if (nlp_opts->print_level > 0)
        printf("ocp_nlp_sqp: performing SOC, since prelim. line search returned %d\n\n", line_search_status);
    int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    // int *nv = dims->nv;
    // int *ni = dims->ni;

    /* evaluate constraints & dynamics at new step */
    // NOTE: setting up the new iterate and evaluating is not needed here,
    //   since this evaluation was perfomed just before this call in the early terminated line search.

    // NOTE: similar to ocp_nlp_evaluate_merit_fun
    // update QP rhs
    // d_i = c_i(x_k + p_k) - \nabla c_i(x_k)^T * p_k
    struct blasfeo_dvec *tmp_fun_vec;

    for (ii = 0; ii <= N; ii++)
    {
        if (ii < N)
        {
            // b -- dynamics
            tmp_fun_vec = config->dynamics[ii]->memory_get_fun_ptr(nlp_mem->dynamics[ii]);
            // add - \nabla c_i(x_k)^T * p_k
            // c_i = f(x_k, u_k) - x_{k+1} (see dynamics module)
            blasfeo_dgemv_t(nx[ii]+nu[ii], nx[ii+1], -1.0, qp_in->BAbt+ii, 0, 0,
                            qp_out->ux+ii, 0, -1.0, tmp_fun_vec, 0, qp_in->b+ii, 0);
            // NOTE: not sure why it is - tmp_fun_vec here!
            blasfeo_dvecad(nx[ii+1], 1.0, qp_out->ux+ii+1, nu[ii+1], qp_in->b+ii, 0);
        }

        /* INEQUALITIES */
        // d -- constraints
        tmp_fun_vec = config->constraints[ii]->memory_get_fun_ptr(nlp_mem->constraints[ii]);
        /* SOC for bounds can be skipped (because linear) */
        // NOTE: SOC can also be skipped for truely linear constraint, i.e. ng of nlp,
        //      now using ng of QP = (nh+ng)

        // upper & lower
        blasfeo_dveccp(ng[ii], tmp_fun_vec, nb[ii], qp_in->d+ii, nb[ii]); // lg
        blasfeo_dveccp(ng[ii], tmp_fun_vec, 2*nb[ii]+ng[ii], qp_in->d+ii, 2*nb[ii]+ng[ii]); // ug
        // general linear / linearized!
        // tmp_ni = D * u + C * x
        blasfeo_dgemv_t(nu[ii]+nx[ii], ng[ii], 1.0, qp_in->DCt+ii, 0, 0, qp_out->ux+ii, 0,
                        0.0, &work->nlp_work->tmp_ni, 0, &work->nlp_work->tmp_ni, 0);
        // d[nb:nb+ng] += tmp_ni (lower)
        blasfeo_dvecad(ng[ii], 1.0, &work->nlp_work->tmp_ni, 0, qp_in->d+ii, nb[ii]);
        // d[nb:nb+ng] -= tmp_ni
        blasfeo_dvecad(ng[ii], -1.0, &work->nlp_work->tmp_ni, 0, qp_in->d+ii, 2*nb[ii]+ng[ii]);

        // add slack contributions
        // d[nb:nb+ng] += slack[idx]
        // qp_in->idxs_rev
        for (int j = 0; j < nb[ii]+ng[ii]; j++)
        {
            int slack_index = qp_in->idxs_rev[ii][j];
            if (slack_index >= 0)
            {
                // add slack contribution for lower and upper constraint
                // lower
                BLASFEO_DVECEL(qp_in->d+ii, j) -=
                        BLASFEO_DVECEL(qp_out->ux+ii, slack_index+nx[ii]+nu[ii]);
                // upper
                BLASFEO_DVECEL(qp_in->d+ii, j+nb[ii]+ng[ii]) -=
                        BLASFEO_DVECEL(qp_out->ux+ii, slack_index+nx[ii]+nu[ii]+ns[ii]);
            }
        }

        // NOTE: bounds on slacks can be skipped, since they are linear.
        // blasfeo_daxpy(2*ns[ii], -1.0, qp_out->ux+ii, nx[ii]+nu[ii], qp_in->d+ii, 2*nb[ii]+2*ng[ii], qp_in->d+ii, 2*nb[ii]+2*ng[ii]);

        // printf("SOC: qp_in->d final value\n");
        // blasfeo_print_exp_dvec(2*nb[ii]+2*ng[ii], qp_in->d+ii, 0);
    }

    if (nlp_opts->print_level > 3)
    {
        printf("\n\nSQP: SOC ocp_qp_in at iteration %d\n", sqp_iter);
        print_ocp_qp_in(qp_in);
    }

#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
    ocp_nlp_sqp_dump_qp_in_to_file(qp_in, sqp_iter, 1);
#endif

    // solve QP
    // acados_tic(&timer1);
    int qp_status = qp_solver->evaluate(qp_solver, dims->qp_solver, qp_in, qp_out,
                                    opts->nlp_opts->qp_solver_opts, nlp_mem->qp_solver_mem, nlp_work->qp_work);
    // NOTE: QP is not timed, since this computation time is attributed to globalization.

    // compute correct dual solution in case of Hessian regularization
    config->regularize->correct_dual_sol(config->regularize, dims->regularize,
                                        opts->nlp_opts->regularize, nlp_mem->regularize_mem);

    ocp_qp_out_get(qp_out, "qp_info", &qp_info_);
    int qp_iter = qp_info_->num_iter;

    // save statistics of last qp solver call
    // TODO: SOC QP solver call should be warm / hot started!
    if (sqp_iter+1 < mem->stat_m)
    {
        // mem->stat[mem->stat_n*(sqp_iter+1)+4] = qp_status;
        // add qp_iter; should maybe be in a seperate statistic
        mem->stat[mem->stat_n*(sqp_iter+1)+5] += qp_iter;
    }

    // compute external QP residuals (for debugging)
    if (opts->ext_qp_res)
    {
        ocp_qp_res_compute(qp_in, qp_out, work->qp_res, work->qp_res_ws);
        if (sqp_iter+1 < mem->stat_m)
            ocp_qp_res_compute_nrm_inf(work->qp_res, mem->stat+(mem->stat_n*(sqp_iter+1)+7));
    }

    if (nlp_opts->print_level > 3)
    {
        printf("\n\nSQP: SOC ocp_qp_out at iteration %d\n", sqp_iter);
        print_ocp_qp_out(qp_out);
    }

#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
        ocp_nlp_sqp_dump_qp_out_to_file(qp_out, sqp_iter, 1);
#endif

    // exit conditions on QP status
    if ((qp_status!=ACADOS_SUCCESS) & (qp_status!=ACADOS_MAXITER))
    {
#ifndef ACADOS_SILENT
        printf("\nQP solver returned error status %d in SQP iteration %d for SOC QP in QP iteration %d.\n",
            qp_status, sqp_iter, qp_iter);
#endif
        if (nlp_opts->print_level > 1)
        {
            printf("\nFailed to solve the following QP:\n");
            if (nlp_opts->print_level > 3)
                print_ocp_qp_in(qp_in);
        }

        mem->status = ACADOS_QP_FAILURE;
        mem->sqp_iter = sqp_iter;

        return ACADOS_QP_FAILURE;
    }
    return true;
}

static int ocp_nlp_line_search_merit_check_full_step(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work, int sqp_iter)
{
    int N = dims->N;
    int *nv = dims->nv;

    ocp_qp_out *qp_out = mem->qp_out;

    double merit_fun1;

    // copy out (current iterate) to work->tmp_nlp_out
    for (int i = 0; i <= N; i++)
        blasfeo_dveccp(nv[i], out->ux+i, 0, work->tmp_nlp_out->ux+i, 0);
    // NOTE: copying duals not needed, as they dont enter the merit function, see ocp_nlp_line_search


    /* modify/initialize merit function weights (Leineweber1999 M5.1, p.89) */
    if (sqp_iter==0)
    {
        merit_backtracking_initialize_weights(dims, work->weight_merit_fun, qp_out);
    }
    else
    {
        // backup weights
        copy_multipliers_nlp_to_qp(dims, work->weight_merit_fun, work->tmp_qp_out);
        // update weights
        merit_backtracking_update_weights(dims, work->weight_merit_fun, qp_out);
    }

    double merit_fun0 = ocp_nlp_evaluate_merit_fun(config, dims, in, out, opts, mem, work);
    double alpha = 1.0;

    // TODO(oj): should the merit weight update be undone in case of early termination?
    double violation_current = ocp_nlp_get_violation_inf_norm(config, dims, in, out, opts, mem, work);

    // tmp_nlp_out = out + alpha * qp_out
    for (int i = 0; i <= N; i++)
        blasfeo_daxpy(nv[i], alpha, qp_out->ux+i, 0, out->ux+i, 0, work->tmp_nlp_out->ux+i, 0);
    merit_fun1 = ocp_nlp_evaluate_merit_fun(config, dims, in, out, opts, mem, work);

    double violation_step = ocp_nlp_get_violation_inf_norm(config, dims, in, out, opts, mem, work);
    if (opts->print_level > 0)
    {
        printf("\npreliminary line_search: merit0 %e, merit1 %e; viol_current %e, viol_step %e\n", merit_fun0, merit_fun1, violation_current, violation_step);
    }

    if (isnan(merit_fun1) || isinf(merit_fun1))
    {
        // do nothing and continue with normal line search, i.e. step reduction
        if (sqp_iter != 0)
        {
            // reset merit function weights;
            copy_multipliers_qp_to_nlp(dims, work->tmp_qp_out, work->weight_merit_fun);
        }
        return ACADOS_NAN_DETECTED;
    }
    if (merit_fun1 < merit_fun0 && violation_step < violation_current)
    {
        // full step if merit and constraint violation improves
        // TODO: check armijo in this case?
        return ACADOS_SUCCESS;
    }
    else
    {
        // trigger SOC
        if (sqp_iter != 0)
        {
            // reset merit function weights;
            copy_multipliers_qp_to_nlp(dims, work->tmp_qp_out, work->weight_merit_fun);
        }
        return ACADOS_MINSTEP;
    }
}

static double ocp_nlp_get_violation_inf_norm(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                  ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
                                  ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    // computes constraint violation infinity norm
    // assumes constraint functions are evaluated before, e.g. done in ocp_nlp_evaluate_merit_fun
    int i, j;
    int N = dims->N;
    int *nx = dims->nx;
    int *ni = dims->ni;
    struct blasfeo_dvec *tmp_fun_vec;
    double violation = 0.0;
    double tmp;
    for (i=0; i<N; i++)
    {
        tmp_fun_vec = config->dynamics[i]->memory_get_fun_ptr(mem->dynamics[i]);
        for (j=0; j<nx[i+1]; j++)
        {
            tmp = fabs(BLASFEO_DVECEL(tmp_fun_vec, j));
            violation = tmp > violation ? tmp : violation;
        }
    }

    for (i=0; i<=N; i++)
    {
        tmp_fun_vec = config->constraints[i]->memory_get_fun_ptr(mem->constraints[i]);
        for (j=0; j<2*ni[i]; j++)
        {
            // Note constraint violation corresponds to > 0
            tmp = BLASFEO_DVECEL(tmp_fun_vec, j);
            violation = tmp > violation ? tmp : violation;
        }
    }

    return violation;
}


static void copy_multipliers_nlp_to_qp(ocp_nlp_dims *dims, ocp_nlp_out *from, ocp_qp_out *to)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *ni = dims->ni;
    for (int i = 0; i <= N; i++)
    {
        blasfeo_dveccp(2*ni[i], from->lam+i, 0, to->lam+i, 0);
    }
    for (int i = 0; i < N; i++)
    {
        blasfeo_dveccp(nx[i+1], from->pi+i, 0, to->pi+i, 0);
    }
    return;
}


static void copy_multipliers_qp_to_nlp(ocp_nlp_dims *dims, ocp_qp_out *from, ocp_nlp_out *to)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *ni = dims->ni;
    for (int i = 0; i <= N; i++)
    {
        blasfeo_dveccp(2*ni[i], from->lam+i, 0, to->lam+i, 0);
    }
    for (int i = 0; i < N; i++)
    {
        blasfeo_dveccp(nx[i+1], from->pi+i, 0, to->pi+i, 0);
    }
    return;
}

int ocp_nlp_globalization_merit_backtracking_find_acceptable_iterate(void *nlp_config_, void *nlp_dims_, void *nlp_in_, void *nlp_out_, void *nlp_mem_, void *nlp_work_, void *nlp_opts_)
{
    ocp_nlp_config *nlp_config = nlp_config_;
    ocp_nlp_dims *nlp_dims = nlp_dims_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = nlp_mem_;
    ocp_nlp_workspace *nlp_work = nlp_work_;
    ocp_nlp_opts *nlp_opts = nlp_opts_;
    
    printf("Merit backtracking line search\n");

    // int sqp_iter = 1;
    // bool do_line_search = true;
//     if (nlp_opts->globalization->globalization_use_SOC)
//     {
//         do_line_search = ocp_nlp_soc_line_search(nlp_config, nlp_dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, sqp_iter);
// //         if (nlp_mem->status == ACADOS_QP_FAILURE)
// //         {
// // #if defined(ACADOS_WITH_OPENMP)
// //             // restore number of threads
// //             omp_set_num_threads(num_threads_bkp);
// // #endif
// //             // mem->time_tot = acados_toc(&timer0);
// //             return mem->status;
// //         }
//     }

    // if (do_line_search)
    // {
    //     int line_search_status;
    //     line_search_status = ocp_nlp_line_search(nlp_config, nlp_dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, sqp_iter, &mem->alpha);
    //     if (line_search_status == ACADOS_NAN_DETECTED)
    //     {
    //         mem->status = ACADOS_NAN_DETECTED;
    //         return mem->status;
    //     }
    // }
    // // mem->time_glob += acados_toc(&timer1);
    // // nlp_mem->stat[mem->stat_n*(sqp_iter+1)+6] = mem->alpha;

    // // update variables
    // ocp_nlp_update_variables_sqp(nlp_config, nlp_dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, nlp_out, mem->alpha);
}

int ocp_nlp_globalization_merit_backtracking_needs_objective_value()
{
    return 1;
}

void ocp_nlp_globalization_merit_backtracking_config_initialize_default(ocp_nlp_globalization_config *config)
{
    // opts
    config->opts_calculate_size = &ocp_nlp_globalization_merit_backtracking_opts_calculate_size;
    config->opts_assign = &ocp_nlp_globalization_merit_backtracking_opts_assign;
    config->opts_initialize_default = &ocp_nlp_globalization_merit_backtracking_opts_initialize_default;
    config->opts_set = &ocp_nlp_globalization_merit_backtracking_opts_set;
    // functions
    config->find_acceptable_iterate = &ocp_nlp_globalization_merit_backtracking_find_acceptable_iterate;
    config->print_iteration_header = &ocp_nlp_globalization_merit_backtracking_print_iteration_header;
    config->print_iteration = &ocp_nlp_globalization_merit_backtracking_print_iteration;
    config->needs_objective_value = &ocp_nlp_globalization_merit_backtracking_needs_objective_value;
}