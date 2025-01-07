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
#include "acados/utils/mem.h"

// blasfeo
#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"


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
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_globalization_merit_backtracking_opts *opts = opts_;
    ocp_nlp_globalization_opts *globalization_opts = opts->globalization_opts;
    ocp_nlp_globalization_config *config = config_;

    ocp_nlp_globalization_opts_initialize_default(config, dims, globalization_opts);
    return;
}


void ocp_nlp_globalization_merit_backtracking_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_globalization_merit_backtracking_opts *opts = opts_;
    ocp_nlp_globalization_config *config = config_;

    ocp_nlp_globalization_opts_set(config, opts->globalization_opts, field, value);

    return;
}

void *ocp_nlp_globalization_merit_backtracking_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_globalization_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_globalization_merit_backtracking_opts *opts = (ocp_nlp_globalization_merit_backtracking_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_globalization_merit_backtracking_opts);

    opts->globalization_opts = ocp_nlp_globalization_opts_assign(config, dims, c_ptr);
    c_ptr += ocp_nlp_globalization_opts_calculate_size(config, dims);

    assert((char *) raw_memory + ocp_nlp_globalization_merit_backtracking_opts_calculate_size(config_, dims_) >=
           c_ptr);

    return opts;
}

/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_globalization_merit_backtracking_memory_calculate_size(void *config_, void *dims_)
{
    acados_size_t size = 0;

    size += sizeof(ocp_nlp_globalization_merit_backtracking_memory);

    return size;
}

void *ocp_nlp_globalization_merit_backtracking_memory_assign(void *config_, void *dims_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // initial align
    align_char_to(8, &c_ptr);

    ocp_nlp_globalization_merit_backtracking_memory *mem = (ocp_nlp_globalization_merit_backtracking_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_globalization_merit_backtracking_memory);

    align_char_to(8, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_globalization_merit_backtracking_memory_calculate_size(config_, dims_) >= c_ptr);

    return mem;
}

/************************************************
 * functions
 ************************************************/

int ocp_nlp_line_search(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work,
            int sqp_iter, double *alpha_reference)
{
    ocp_nlp_globalization_merit_backtracking_opts *merit_opts = opts->globalization;
    ocp_nlp_globalization_opts *globalization_opts = merit_opts->globalization_opts;
    int i, j;

    int N = dims->N;
    int *nv = dims->nv;

    double merit_fun1 = 0;
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

void ocp_nlp_globalization_merit_backtracking_print_iteration(double objective_value,
                                                            int iter_count,
                                                            void* nlp_res_,
                                                            double step_norm,
                                                            double reg_param,
                                                            int qp_status,
                                                            int qp_iter,
                                                            void* nlp_opts_,
                                                            void* mem_)
{
    ocp_nlp_res *nlp_res = nlp_res_;
    ocp_nlp_globalization_merit_backtracking_memory* mem = mem_;

    if ((iter_count % 10 == 0)){
        ocp_nlp_globalization_merit_backtracking_print_iteration_header();
    }
    printf("%i\t%e\t%e\t%e\t%e\t%d\t%d\t%e\n",
        iter_count,
        nlp_res->inf_norm_res_stat,
        nlp_res->inf_norm_res_eq,
        nlp_res->inf_norm_res_ineq,
        nlp_res->inf_norm_res_comp,
        qp_status,
        qp_iter,
        mem->alpha);
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

static bool ocp_nlp_soc_line_search(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in,
            ocp_nlp_out *nlp_out, ocp_nlp_opts *nlp_opts, ocp_nlp_memory *nlp_mem, ocp_nlp_workspace *nlp_work, int sqp_iter)
{
    int ii;
    int N = dims->N;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_qp_in *qp_in = nlp_mem->qp_in;
    ocp_qp_out *qp_out = nlp_mem->qp_out;
    ocp_nlp_globalization_merit_backtracking_memory *merit_mem = nlp_mem->globalization;
    // qp_info *qp_info_;
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
        merit_mem->alpha = 1.0;
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
                        0.0, &nlp_work->tmp_ni, 0, &nlp_work->tmp_ni, 0);
        // d[nb:nb+ng] += tmp_ni (lower)
        blasfeo_dvecad(ng[ii], 1.0, &nlp_work->tmp_ni, 0, qp_in->d+ii, nb[ii]);
        // d[nb:nb+ng] -= tmp_ni
        blasfeo_dvecad(ng[ii], -1.0, &nlp_work->tmp_ni, 0, qp_in->d+ii, 2*nb[ii]+ng[ii]);

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
        // print_ocp_qp_in(qp_in);
    }

#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
    ocp_nlp_dump_qp_in_to_file(qp_in, sqp_iter, 1);
#endif

    // solve QP
    // acados_tic(&timer1);
    int qp_status = qp_solver->evaluate(qp_solver, dims->qp_solver, qp_in, qp_out,
                                    nlp_opts->qp_solver_opts, nlp_mem->qp_solver_mem, nlp_work->qp_work);
    // NOTE: QP is not timed, since this computation time is attributed to globalization.

    // compute correct dual solution in case of Hessian regularization
    config->regularize->correct_dual_sol(config->regularize, dims->regularize,
                                        nlp_opts->regularize, nlp_mem->regularize_mem);

    // ocp_qp_out_get(qp_out, "qp_info", &qp_info_);
    // int qp_iter = qp_info_->num_iter;

    // save statistics of last qp solver call
    // TODO: SOC QP solver call should be warm / hot started!
    // if (sqp_iter+1 < nlp_mem->stat_m)
    // {
    //     // mem->stat[mem->stat_n*(sqp_iter+1)+4] = qp_status;
    //     // add qp_iter; should maybe be in a seperate statistic
    //     nlp_mem->stat[nlp_mem->stat_n*(sqp_iter+1)+5] += qp_iter;
    // }

    // compute external QP residuals (for debugging)
    // if (nlp_opts->ext_qp_res)
    // {
    //     ocp_qp_res_compute(qp_in, qp_out, nlp_work->qp_res, nlp_work->qp_res_ws);
    //     if (sqp_iter+1 < nlp_mem->stat_m)
    //         ocp_qp_res_compute_nrm_inf(nlp_work->qp_res, nlp_mem->stat+(nlp_mem->stat_n*(sqp_iter+1)+7));
    // }

    // if (nlp_opts->print_level > 3)
    // {
    //     printf("\n\nSQP: SOC ocp_qp_out at iteration %d\n", sqp_iter);
    //     print_ocp_qp_out(qp_out);
    // }

#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
        ocp_nlp_dump_qp_out_to_file(qp_out, sqp_iter, 1);
#endif

    // exit conditions on QP status
    if ((qp_status!=ACADOS_SUCCESS) & (qp_status!=ACADOS_MAXITER))
    {
#ifndef ACADOS_SILENT
        printf("\nQP solver returned error status %d in SQP iteration %d for SOC QP.\n",
            qp_status, sqp_iter);
#endif
        // if (nlp_opts->print_level > 1)
        // {
        //     printf("\nFailed to solve the following QP:\n");
        //     if (nlp_opts->print_level > 3)
        //         print_ocp_qp_in(qp_in);
        // }

        nlp_mem->status = ACADOS_QP_FAILURE;
        nlp_mem->iter = sqp_iter;

        return true;
    }
    return true;
}

double ocp_nlp_compute_merit_gradient(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                  ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
                                  ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    /* computes merit function gradient at iterate: out -- using already evaluated gradients of submodules
       with weights: work->weight_merit_fun */
    int i, j;

    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;

    double merit_grad = 0.0;
    double weight;

    // NOTE: step is in: mem->qp_out->ux
    struct blasfeo_dvec *tmp_vec; // size nv
    struct blasfeo_dvec tmp_vec_nxu = work->tmp_nv;  // size nxu
    struct blasfeo_dvec dxnext_dy = work->dxnext_dy;  // size nx

    // cost
    for (i=0; i<=N; i++)
    {
        tmp_vec = config->cost[i]->memory_get_grad_ptr(mem->cost[i]);
        merit_grad += blasfeo_ddot(nv[i], tmp_vec, 0, mem->qp_out->ux + i, 0);
    }
    double merit_grad_cost = merit_grad;

    /* dynamics */
    double merit_grad_dyn = 0.0;
    for (i=0; i<N; i++)
    {
        // get shooting node gap x_next(x_n, u_n) - x_{n+1};
        tmp_vec = config->dynamics[i]->memory_get_fun_ptr(mem->dynamics[i]);

        /* compute directional derivative of xnext with direction y -> dxnext_dy */
        blasfeo_dgemv_t(nx[i]+nu[i], nx[i+1], 1.0, mem->qp_in->BAbt+i, 0, 0, mem->qp_out->ux+i, 0,
                        0.0, &dxnext_dy, 0, &dxnext_dy, 0);

        /* add merit gradient contributions depending on sign of shooting gap */
        for (j = 0; j < nx[i+1]; j++)
        {
            weight = BLASFEO_DVECEL(work->weight_merit_fun->pi+i, j);
            double deqj_dy = BLASFEO_DVECEL(&dxnext_dy, j) - BLASFEO_DVECEL(mem->qp_out->ux+(i+1), nu[i+1]+j);
            {
                if (BLASFEO_DVECEL(tmp_vec, j) > 0)
                {
                    merit_grad_dyn += weight * deqj_dy;
                    // printf("\ndyn_contribution +%e, weight %e, deqj_dy %e, i %d, j %d", weight * deqj_dy, weight, deqj_dy, i, j);
                }
                else
                {
                    merit_grad_dyn -= weight * deqj_dy;
                    // printf("\ndyn_contribution %e, weight %e, deqj_dy %e, i %d, j %d", -weight * deqj_dy, weight, deqj_dy, i, j);
                }
            }
        }
    }

    /* inequality contributions */
    // NOTE: slack bound inequalities are not considered here.
    // They should never be infeasible. Only if explicitly initialized infeasible from outside.
    int constr_index, slack_index_in_ux, slack_index;
    ocp_qp_dims* qp_dims = mem->qp_in->dim;
    int *nb = qp_dims->nb;
    int *ng = qp_dims->ng;
    int *ns = qp_dims->ns;
    double merit_grad_ineq = 0.0;
    double slack_step;

    for (i=0; i<=N; i++)
    {
        tmp_vec = config->constraints[i]->memory_get_fun_ptr(mem->constraints[i]);
        int *idxb = mem->qp_in->idxb[i];
        if (ni[i] > 0)
        {
            // NOTE: loop could be simplified handling lower and upper constraints together.
            for (j = 0; j < 2 * (nb[i] + ng[i]); j++) // 2 * ni
            {
                double constraint_val = BLASFEO_DVECEL(tmp_vec, j);
                if (constraint_val > 0)
                {
                    weight = BLASFEO_DVECEL(work->weight_merit_fun->lam+i, j);

                    // find corresponding slack value
                    constr_index = j < nb[i]+ng[i] ? j : j-(nb[i]+ng[i]);
                    slack_index = mem->qp_in->idxs_rev[i][constr_index];
                    // if softened: add slack contribution
                    if (slack_index >= 0)
                    {
                        slack_index_in_ux = j < (nb[i]+ng[i]) ? nx[i] + nu[i] + slack_index
                                                              : nx[i] + nu[i] + slack_index + ns[i];
                        slack_step = BLASFEO_DVECEL(mem->qp_out->ux+i, slack_index_in_ux);
                        merit_grad_ineq -= weight * slack_step;
                        // printf("at node %d, ineq %d, idxs_rev[%d] = %d\n", i, j, constr_index, slack_index);
                        // printf("slack contribution: uxs[%d] = %e\n", slack_index_in_ux, slack_step);
                    }


                    // NOTE: the inequalities are internally organized in the following order:
                    //     [ lbu lbx lg lh lphi ubu ubx ug uh uphi;
                    //     lsbu lsbx lsg lsh lsphi usbu usbx usg ush usphi]
                    // printf("constraint %d %d is active with value %e", i, j, constraint_val);
                    if (j < nb[i])
                    {
                        // printf("lower idxb[%d] = %d dir %f, constraint_val %f, nb = %d\n", j, idxb[j], BLASFEO_DVECEL(mem->qp_out->ux, idxb[j]), constraint_val, nb[i]);
                        merit_grad_ineq += weight * BLASFEO_DVECEL(mem->qp_out->ux+i, idxb[j]);
                    }
                    else if (j < nb[i] + ng[i])
                    {
                        // merit_grad_ineq += weight * mem->qp_in->DCt_j * dux
                        blasfeo_dcolex(nx[i] + nu[i], mem->qp_in->DCt+i, j - nb[i], 0, &tmp_vec_nxu, 0);
                        merit_grad_ineq += weight * blasfeo_ddot(nx[i] + nu[i], &tmp_vec_nxu, 0, mem->qp_out->ux+i, 0);
                        // printf("general linear constraint lower contribution = %e, val = %e\n", blasfeo_ddot(nx[i] + nu[i], &tmp_vec_nxu, 0, mem->qp_out->ux+i, 0), constraint_val);
                    }
                    else if (j < 2*nb[i] + ng[i])
                    {
                        // printf("upper idxb[%d] = %d dir %f, constraint_val %f, nb = %d\n", j-nb[i]-ng[i], idxb[j-nb[i]-ng[i]], BLASFEO_DVECEL(mem->qp_out->ux, idxb[j-nb[i]-ng[i]]), constraint_val, nb[i]);
                        merit_grad_ineq += weight * BLASFEO_DVECEL(mem->qp_out->ux+i, idxb[j-nb[i]-ng[i]]);
                    }
                    else if (j < 2*nb[i] + 2*ng[i])
                    {
                        blasfeo_dcolex(nx[i] + nu[i], mem->qp_in->DCt+i, j - 2*nb[i] - ng[i], 0, &tmp_vec_nxu, 0);
                        merit_grad_ineq += weight * blasfeo_ddot(nx[i] + nu[i], &tmp_vec_nxu, 0, mem->qp_out->ux+i, 0);
                        // printf("general linear constraint upper contribution = %e, val = %e\n", blasfeo_ddot(nx[i] + nu[i], &tmp_vec_nxu, 0, mem->qp_out->ux+i, 0), constraint_val);
                    }
                }
            }
        }
    }
    // print_ocp_qp_dims(qp_dims);
    // print_ocp_qp_in(mem->qp_in);

    merit_grad = merit_grad_cost + merit_grad_dyn + merit_grad_ineq;
    if (opts->print_level > 1)
        printf("computed merit_grad = %e, merit_grad_cost = %e, merit_grad_dyn = %e, merit_grad_ineq = %e\n", merit_grad, merit_grad_cost, merit_grad_dyn, merit_grad_ineq);

    return merit_grad;
}



double ocp_nlp_evaluate_merit_fun(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                  ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
                                  ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    /* computes merit function value at iterate: tmp_nlp_out, with weights: work->weight_merit_fun */
    //int j;

    int N = dims->N;
    int *nx = dims->nx;
    int *ni = dims->ni;

    double merit_fun = 0.0;

    // set evaluation point to tmp_nlp_out
    ocp_nlp_set_primal_variable_pointers_in_submodules(config, dims, in, work->tmp_nlp_out, mem);
    // compute fun value
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i=0; i<N; i++)
    {
        // dynamics: Note has to be first, because cost_integration might be used.
        config->dynamics[i]->compute_fun(config->dynamics[i], dims->dynamics[i], in->dynamics[i],
                                         opts->dynamics[i], mem->dynamics[i], work->dynamics[i]);
    }
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i=0; i<=N; i++)
    {
        // cost
        config->cost[i]->compute_fun(config->cost[i], dims->cost[i], in->cost[i], opts->cost[i],
                                    mem->cost[i], work->cost[i]);
    }
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i=0; i<=N; i++)
    {
        // constr
        config->constraints[i]->compute_fun(config->constraints[i], dims->constraints[i],
                                            in->constraints[i], opts->constraints[i],
                                            mem->constraints[i], work->constraints[i]);
    }
    // reset evaluation point to SQP iterate
    ocp_nlp_set_primal_variable_pointers_in_submodules(config, dims, in, out, mem);

    double *tmp_fun;
    double tmp;
    struct blasfeo_dvec *tmp_fun_vec;

    double cost_fun = 0.0;
    for(int i=0; i<=N; i++)
    {
        tmp_fun = config->cost[i]->memory_get_fun_ptr(mem->cost[i]);
        cost_fun += *tmp_fun;
    }

    double dyn_fun = 0.0;
    for(int i=0; i<N; i++)
    {
        tmp_fun_vec = config->dynamics[i]->memory_get_fun_ptr(mem->dynamics[i]);
        // printf("\nMerit: dyn will multiply tmp_fun, weights %d\n", i);
        // blasfeo_print_exp_tran_dvec(nx[i+1], tmp_fun_vec, 0);
        // blasfeo_print_exp_tran_dvec(nx[i+1], work->weight_merit_fun->pi+i, 0);
        for(int j=0; j<nx[i+1]; j++)
        {
//            printf("\n%e %e\n", fabs(BLASFEO_DVECEL(work->weight_merit_fun->pi+i, j)), fabs(BLASFEO_DVECEL(tmp_fun_vec, j)));
            dyn_fun += fabs(BLASFEO_DVECEL(work->weight_merit_fun->pi+i, j)) * fabs(BLASFEO_DVECEL(tmp_fun_vec, j));
        }
    }

    double constr_fun = 0.0;
    for(int i=0; i<=N; i++)
    {
//        printf("\ni %d\n", i);
        tmp_fun_vec = config->constraints[i]->memory_get_fun_ptr(mem->constraints[i]);
//        blasfeo_print_exp_tran_dvec(2*ni[i], tmp_fun_vec, 0);
//        blasfeo_print_exp_tran_dvec(2*ni[i], work->weight_merit_fun->lam+i, 0);
        for (int j=0; j<2*ni[i]; j++)
        {
            tmp = BLASFEO_DVECEL(tmp_fun_vec, j);
            if (tmp > 0.0)
            {
                // tmp = constraint violation
                // printf("IN merit fun: ineq i %d, j %d tmp_fun %e, multiplier %e\n", i, j, tmp, BLASFEO_DVECEL(work->weight_merit_fun->lam+i, j));
                constr_fun += fabs(BLASFEO_DVECEL(work->weight_merit_fun->lam+i, j)) * tmp;
            }
        }
    }

    merit_fun = cost_fun + dyn_fun + constr_fun;

    // printf("Merit fun: %e cost: %e dyn: %e constr: %e\n", merit_fun, cost_fun, dyn_fun, constr_fun);

    return merit_fun;
}


void merit_backtracking_initialize_weights(ocp_nlp_dims *dims, ocp_nlp_out *weight_merit_fun, ocp_qp_out *qp_out)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *ni = dims->ni;
    // equality merit weights = abs( eq multipliers of qp_sol )
    for (int i = 0; i < N; i++)
    {
        for (int j=0; j<nx[i+1]; j++)
        {
            BLASFEO_DVECEL(weight_merit_fun->pi+i, j) = fabs(BLASFEO_DVECEL(qp_out->pi+i, j));
        }
    }

    for (int i = 0; i <= N; i++)
    {
        blasfeo_dveccp(2*ni[i], qp_out->lam+i, 0, weight_merit_fun->lam+i, 0);
    }
}

void merit_backtracking_update_weights(ocp_nlp_dims *dims, ocp_nlp_out *weight_merit_fun, ocp_qp_out *qp_out)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *ni = dims->ni;
    double tmp0, tmp1;

    // update weights
    for (int i = 0; i < N; i++)
    {
        for (int j=0; j<nx[i+1]; j++)
        {
            // abs(lambda) (LW)
            tmp0 = fabs(BLASFEO_DVECEL(qp_out->pi+i, j));
            // .5 * (abs(lambda) + sigma)
            tmp1 = 0.5 * (tmp0 + BLASFEO_DVECEL(weight_merit_fun->pi+i, j));
            BLASFEO_DVECEL(weight_merit_fun->pi+i, j) = tmp0 > tmp1 ? tmp0 : tmp1;
        }
    }
    for (int i = 0; i <= N; i++)
    {
        for (int j=0; j<2*ni[i]; j++)
        {
            // mu (LW)
            tmp0 = BLASFEO_DVECEL(qp_out->lam+i, j);
            // .5 * (mu + tau)
            tmp1 = 0.5 * (tmp0 + BLASFEO_DVECEL(weight_merit_fun->lam+i, j));
            BLASFEO_DVECEL(weight_merit_fun->lam+i, j) = tmp0>tmp1 ? tmp0 : tmp1;
        }
    }
}


static int ocp_nlp_ddp_backtracking_line_search(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out,
                ocp_nlp_memory *nlp_mem, void* solver_mem, ocp_nlp_workspace *nlp_work, ocp_nlp_opts *nlp_opts)
{
    // evaluate the objective of the QP (as predicted reduction)
    ocp_nlp_globalization_merit_backtracking_opts *merit_opts = nlp_opts->globalization;
    ocp_nlp_globalization_opts *globalization_opts = merit_opts->globalization_opts;
    ocp_nlp_globalization_merit_backtracking_memory *mem = nlp_mem->globalization;
    int N = dims->N;
    double pred = -nlp_mem->qp_cost_value;
    double alpha = 1.0;
    double trial_cost;
    double negative_ared;
    double *tmp_fun;

    int i;

    while (true)
    {
        // Do the DDP forward sweep to get the trial iterate
        config->step_update(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem,
                                     nlp_work, nlp_work->tmp_nlp_out, solver_mem, alpha, globalization_opts->full_step_dual);

        ///////////////////////////////////////////////////////////////////////
        // Evaluate cost function at trial iterate
        // set evaluation point to tmp_nlp_out
        ocp_nlp_set_primal_variable_pointers_in_submodules(config, dims, nlp_in, nlp_work->tmp_nlp_out, nlp_mem);
        // compute fun value
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
        for (i=0; i<=N; i++)
        {
            // cost
            config->cost[i]->compute_fun(config->cost[i], dims->cost[i], nlp_in->cost[i], nlp_opts->cost[i],
                                        nlp_mem->cost[i], nlp_work->cost[i]);
        }
        ocp_nlp_set_primal_variable_pointers_in_submodules(config, dims, nlp_in, nlp_out, nlp_mem);
        trial_cost = 0.0;
        for(i=0; i<=N; i++)
        {
            tmp_fun = config->cost[i]->memory_get_fun_ptr(nlp_mem->cost[i]);
            trial_cost += *tmp_fun;
        }

        negative_ared = trial_cost - nlp_mem->cost_value;
        // Check Armijo sufficient decrease condition
        if (negative_ared <= fmin(-globalization_opts->eps_sufficient_descent*alpha* fmax(pred, 0) + 1e-18, 0))
        {
            // IF step accepted: update x
            // reset evaluation point to SQP iterate
            mem->alpha = alpha;
            nlp_mem->cost_value = trial_cost;
            return ACADOS_SUCCESS;
        }
        else
        {
            // Reduce step size
            alpha *= globalization_opts->alpha_reduction;
        }

        if (alpha < globalization_opts->alpha_min)
        {
            printf("Linesearch: Step size gets too small. Increasing regularization.\n");
            mem->alpha = 0.0; // set to zero such that regularization is increased
            return ACADOS_MINSTEP;
        }
    }
}

int ocp_nlp_globalization_merit_backtracking_find_acceptable_iterate_for_ddp(void *nlp_config_, void *nlp_dims_, void *nlp_in_, void *nlp_out_, void *nlp_mem_, void *solver_mem, void *nlp_work_, void *nlp_opts_, double *step_size)
{
    ocp_nlp_config *nlp_config = nlp_config_;
    ocp_nlp_dims *nlp_dims = nlp_dims_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = nlp_mem_;
    ocp_nlp_globalization_merit_backtracking_memory *mem = nlp_mem->globalization;
    ocp_nlp_workspace *nlp_work = nlp_work_;
    ocp_nlp_opts *nlp_opts = nlp_opts_;

    int linesearch_success = 1;
    // Do the globalization here: Either fixed step or Armijo line search
    // NOTE on timings: currently all within globalization is accounted for within time_glob.
    //   QP solver times could be also attributed there alternatively. Cleanest would be to save them seperately.
    // do backtracking line search on objective function
    linesearch_success = ocp_nlp_ddp_backtracking_line_search(nlp_config, nlp_dims, nlp_in, nlp_out,
                nlp_mem, solver_mem, nlp_work, nlp_opts);

    // Copy new iterate to nlp_out
    if (linesearch_success == ACADOS_SUCCESS)
    {
        // in case line search fails, we do not want to copy trial iterates!
        *step_size = mem->alpha;
        copy_ocp_nlp_out(nlp_dims, nlp_work->tmp_nlp_out, nlp_out);
        return ACADOS_SUCCESS;
    }
    return ACADOS_MINSTEP;
}


int ocp_nlp_globalization_merit_backtracking_find_acceptable_iterate(void *nlp_config_, void *nlp_dims_, void *nlp_in_, void *nlp_out_, void *nlp_mem_, void *solver_mem, void *nlp_work_, void *nlp_opts_, double *step_size)
{
    ocp_nlp_config *nlp_config = nlp_config_;
    ocp_nlp_dims *nlp_dims = nlp_dims_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = nlp_mem_;
    ocp_nlp_globalization_merit_backtracking_memory *mem = nlp_mem->globalization;
    ocp_nlp_workspace *nlp_work = nlp_work_;
    ocp_nlp_opts *nlp_opts = nlp_opts_;
    ocp_nlp_globalization_merit_backtracking_opts *merit_opts = nlp_opts->globalization;
    ocp_nlp_globalization_opts *globalization_opts = merit_opts->globalization_opts;

    int sqp_iter = 1; // NEEDS TO BE CHANGED HERE
    bool do_line_search = true;

    if (merit_opts->globalization_opts->use_SOC)
    {
        do_line_search = ocp_nlp_soc_line_search(nlp_config, nlp_dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, sqp_iter);
        if (nlp_mem->status == ACADOS_QP_FAILURE)
        {
            return nlp_mem->status;
        }
    }

    if (do_line_search)
    {
        int line_search_status;
        line_search_status = ocp_nlp_line_search(nlp_config, nlp_dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, sqp_iter, &mem->alpha);
        if (line_search_status == ACADOS_NAN_DETECTED)
        {
            nlp_mem->status = ACADOS_NAN_DETECTED;
            return nlp_mem->status;
        }
    }

    // update variables
    nlp_config->step_update(nlp_config, nlp_dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, nlp_out, solver_mem, mem->alpha, globalization_opts->full_step_dual);
    *step_size = mem->alpha;
    return ACADOS_SUCCESS;
}

int ocp_nlp_globalization_merit_backtracking_needs_objective_value()
{
    return 1;
}

int ocp_nlp_globalization_merit_backtracking_needs_qp_objective_value()
{
    return 0;
}

int ocp_nlp_globalization_merit_backtracking_ddp_needs_qp_objective_value()
{
    return 1;
}

void ocp_nlp_globalization_merit_backtracking_initialize_memory(void *config_,
    void *dims_, void *nlp_mem_, void *nlp_opts_)
{
    return;
}

void ocp_nlp_globalization_merit_backtracking_config_initialize_default(ocp_nlp_globalization_config *config)
{
    // opts
    config->opts_calculate_size = &ocp_nlp_globalization_merit_backtracking_opts_calculate_size;
    config->opts_assign = &ocp_nlp_globalization_merit_backtracking_opts_assign;
    config->opts_initialize_default = &ocp_nlp_globalization_merit_backtracking_opts_initialize_default;
    config->opts_set = &ocp_nlp_globalization_merit_backtracking_opts_set;
    // memory
    config->memory_calculate_size = &ocp_nlp_globalization_merit_backtracking_memory_calculate_size;
    config->memory_assign = &ocp_nlp_globalization_merit_backtracking_memory_assign;

    // functions
    config->find_acceptable_iterate = &ocp_nlp_globalization_merit_backtracking_find_acceptable_iterate;
    config->print_iteration_header = &ocp_nlp_globalization_merit_backtracking_print_iteration_header;
    config->print_iteration = &ocp_nlp_globalization_merit_backtracking_print_iteration;
    config->needs_objective_value = &ocp_nlp_globalization_merit_backtracking_needs_objective_value;
    config->needs_qp_objective_value = &ocp_nlp_globalization_merit_backtracking_needs_qp_objective_value;
    config->initialize_memory = &ocp_nlp_globalization_merit_backtracking_initialize_memory;
}
