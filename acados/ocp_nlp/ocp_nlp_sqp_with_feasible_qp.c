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


#include "acados/ocp_nlp/ocp_nlp_sqp_with_feasible_qp.h"

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
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_reg_common.h"
#include "acados/ocp_nlp/ocp_nlp_globalization_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados_c/ocp_qp_interface.h"



/************************************************
 * options
 ************************************************/

acados_size_t ocp_nlp_sqp_wfqp_opts_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_sqp_wfqp_opts);

    size += ocp_nlp_opts_calculate_size(config, dims);

    return size;
}



void *ocp_nlp_sqp_wfqp_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_sqp_wfqp_opts *opts = (ocp_nlp_sqp_wfqp_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_sqp_wfqp_opts);

    opts->nlp_opts = ocp_nlp_opts_assign(config, dims, c_ptr);
    c_ptr += ocp_nlp_opts_calculate_size(config, dims);

    assert((char *) raw_memory + ocp_nlp_sqp_wfqp_opts_calculate_size(config, dims) >= c_ptr);

    return opts;
}



void ocp_nlp_sqp_wfqp_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;

    // this first !!!
    ocp_nlp_opts_initialize_default(config, dims, nlp_opts);

    // SQP opts
    opts->max_iter = 20;
    opts->tol_stat = 1e-8;
    opts->tol_eq   = 1e-8;
    opts->tol_ineq = 1e-8;
    opts->tol_comp = 1e-8;
    opts->tol_unbounded = -1e10;
    opts->tol_min_step_norm = 1e-12;

    opts->qp_warm_start = 0;
    opts->warm_start_first_qp = false;
    opts->eval_residual_at_max_iter = false;

    // overwrite default submodules opts
    // qp tolerance
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_stat", &opts->tol_stat);
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_eq", &opts->tol_eq);
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_ineq", &opts->tol_ineq);
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_comp", &opts->tol_comp);

    return;
}



void ocp_nlp_sqp_wfqp_opts_update(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_nlp_opts_update(config, dims, nlp_opts);

    return;
}



void ocp_nlp_sqp_wfqp_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_opts *opts = (ocp_nlp_sqp_wfqp_opts *) opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    int ii;

    char module[MAX_STR_LEN];
    char *ptr_module = NULL;
    int module_length = 0;

    // extract module name
    char *char_ = strchr(field, '_');
    if (char_!=NULL)
    {
        module_length = char_-field;
        for (ii=0; ii<module_length; ii++)
            module[ii] = field[ii];
        module[module_length] = '\0'; // add end of string
        ptr_module = module;
    }

    // pass options to QP module
    if ( ptr_module!=NULL && (!strcmp(ptr_module, "qp")) )
    {
        ocp_nlp_opts_set(config, nlp_opts, field, value);

        if (!strcmp(field, "qp_warm_start"))
        {
            int* i_ptr = (int *) value;
            opts->qp_warm_start = *i_ptr;
        }
    }
    else // nlp opts
    {
        if (!strcmp(field, "max_iter"))
        {
            int* max_iter = (int *) value;
            opts->max_iter = *max_iter;
        }
        else if (!strcmp(field, "tol_stat"))
        {
            double* tol_stat = (double *) value;
            opts->tol_stat = *tol_stat;
            // TODO: set accuracy of the qp_solver to the minimum of current QP accuracy and the one specified.
            config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts, "tol_stat", value);
        }
        else if (!strcmp(field, "tol_eq"))
        {
            double* tol_eq = (double *) value;
            opts->tol_eq = *tol_eq;
            // TODO: set accuracy of the qp_solver to the minimum of current QP accuracy and the one specified.
            config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts, "tol_eq", value);
        }
        else if (!strcmp(field, "tol_ineq"))
        {
            double* tol_ineq = (double *) value;
            opts->tol_ineq = *tol_ineq;
            // TODO: set accuracy of the qp_solver to the minimum of current QP accuracy and the one specified.
            config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts, "tol_ineq", value);
        }
        else if (!strcmp(field, "tol_comp"))
        {
            double* tol_comp = (double *) value;
            opts->tol_comp = *tol_comp;
            // TODO: set accuracy of the qp_solver to the minimum of current QP accuracy and the one specified.
            config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts, "tol_comp", value);
        }
        else if (!strcmp(field, "tol_min_step_norm"))
        {
            double* tol_min_step_norm = (double *) value;
            opts->tol_min_step_norm = *tol_min_step_norm;
        }
        else if (!strcmp(field, "warm_start_first_qp"))
        {
            bool* warm_start_first_qp = (bool *) value;
            opts->warm_start_first_qp = *warm_start_first_qp;
        }
        else if (!strcmp(field, "eval_residual_at_max_iter"))
        {
            bool* eval_residual_at_max_iter = (bool *) value;
            opts->eval_residual_at_max_iter = *eval_residual_at_max_iter;
        }
        else
        {
            ocp_nlp_opts_set(config, nlp_opts, field, value);
        }
    }
    return;
}

void ocp_nlp_sqp_wfqp_opts_set_at_stage(void *config_, void *opts_, size_t stage, const char *field, void* value)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_opts *opts = (ocp_nlp_sqp_wfqp_opts *) opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_nlp_opts_set_at_stage(config, nlp_opts, stage, field, value);

    return;

}

void ocp_nlp_sqp_wfqp_opts_get(void *config_, void *dims_, void *opts_,
                          const char *field, void *return_value_)
{
    // ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_opts *opts = opts_;

    if (!strcmp("nlp_opts", field))
    {
        void **value = return_value_;
        *value = opts->nlp_opts;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_sqp_wfqp_opts_get\n", field);
        exit(1);
    }
}

/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_sqp_wfqp_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_sqp_wfqp_memory);

    // nlp mem
    size += ocp_nlp_memory_calculate_size(config, dims, nlp_opts);

    // primal step norm
    if (opts->nlp_opts->log_primal_step_norm)
    {
        size += opts->max_iter*sizeof(double);
    }
    // stat
    int stat_m = opts->max_iter+1;
    int stat_n = 7;
    if (nlp_opts->ext_qp_res)
        stat_n += 4;
    size += stat_n*stat_m*sizeof(double);

    // idxns
    size += (dims->N + 1) * sizeof(int *);
    int nns, nsbu, nbu, nsbx, nbx, n_nominal_ineq_nlp;
    for (int stage = 0; stage <= dims->N; stage++)
    {
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage], "nsbu", &nsbu);
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage], "nbu", &nbu);
        n_nominal_ineq_nlp = dims->ni[stage] - dims->ns[stage];
        if (stage == 0)
        {
            config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage], "nsbx", &nsbx);
            config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage], "nbx", &nbx);
            nns = n_nominal_ineq_nlp - dims->ns[stage] - nbu + nsbu - nbx + nsbx;
        }
        else
        {
            // We want number of constraints minus number slacked constraints, i.e, ni - ns
            // Then we do not want to slack the control bounds, and we want to avoid the slacked bounds
            nns = n_nominal_ineq_nlp - dims->ns[stage] - nbu + nsbu;
        }
        size += nns * sizeof(int);
    }

    size += 3*8;  // align

    make_int_multiple_of(8, &size);

    return size;
}

void *ocp_nlp_sqp_wfqp_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    char *c_ptr = (char *) raw_memory;

    int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;
    // int *nz = dims->nz;

    // initial align
    align_char_to(8, &c_ptr);

    ocp_nlp_sqp_wfqp_memory *mem = (ocp_nlp_sqp_wfqp_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_sqp_wfqp_memory);

    // assign pointers

    align_char_to(8, &c_ptr);

    // nlp mem
    mem->nlp_mem = ocp_nlp_memory_assign(config, dims, nlp_opts, c_ptr);
    c_ptr += ocp_nlp_memory_calculate_size(config, dims, nlp_opts);

    // primal step norm
    if (opts->nlp_opts->log_primal_step_norm)
    {
        mem->primal_step_norm = (double *) c_ptr;
        c_ptr += opts->max_iter*sizeof(double);
    }

    // stat
    mem->stat_m = opts->max_iter+1;
    mem->stat_n = 7;
    if (nlp_opts->ext_qp_res)
        mem->stat_n += 4;
    mem->stat = (double *) c_ptr;
    assign_and_advance_double(mem->stat_m*mem->stat_n, &mem->stat, &c_ptr);

    align_char_to(8, &c_ptr);

    // ptrs
    assign_and_advance_int_ptrs(N+1, &mem->idxns, &c_ptr);

    // integers
    assign_and_advance_int(N+1, &mem->nns, &c_ptr);

    int nns, nsbu, nbu, nsbx, nbx, n_nominal_ineq_nlp;
    for (int stage = 0; stage <= dims->N; stage++)
    {
        n_nominal_ineq_nlp = dims->ni[stage] - dims->ns[stage];
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage], "nsbu", &nsbu);
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage], "nbu", &nbu);
        if (stage == 0)
        {
            config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage], "nsbx", &nsbx);
            config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage], "nbx", &nbx);
            mem->nns[stage] = n_nominal_ineq_nlp - dims->ns[stage] - nbu + nsbu - nbx + nsbx;
        }
        else
        {
            mem->nns[stage] = n_nominal_ineq_nlp - dims->ns[stage] - nbu + nsbu;
        }
        assign_and_advance_int(mem->nns[stage], &mem->idxns[stage], &c_ptr);
    }

    mem->nlp_mem->status = ACADOS_READY;

    align_char_to(8, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_sqp_wfqp_memory_calculate_size(config, dims, opts) >= c_ptr);

    return mem;
}

/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_nlp_sqp_wfqp_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    acados_size_t size = 0;

    // sqp
    size += sizeof(ocp_nlp_sqp_wfqp_workspace);

    // nlp
    size += ocp_nlp_workspace_calculate_size(config, dims, nlp_opts);

    if (nlp_opts->ext_qp_res)
    {
        // qp res
        size += ocp_qp_res_calculate_size(dims->qp_solver->orig_dims);

        // qp res ws
        size += ocp_qp_res_workspace_calculate_size(dims->qp_solver->orig_dims);
    }

    return size;
}



static void ocp_nlp_sqp_wfqp_cast_workspace(ocp_nlp_config *config, ocp_nlp_dims *dims,
         ocp_nlp_sqp_wfqp_opts *opts, ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_sqp_wfqp_workspace *work)
{
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    // sqp
    char *c_ptr = (char *) work;
    c_ptr += sizeof(ocp_nlp_sqp_wfqp_workspace);

    // nlp
    work->nlp_work = ocp_nlp_workspace_assign(config, dims, nlp_opts, nlp_mem, c_ptr);
    c_ptr += ocp_nlp_workspace_calculate_size(config, dims, nlp_opts);

    if (nlp_opts->ext_qp_res)
    {
        // qp res
        work->nlp_work->qp_res = ocp_qp_res_assign(dims->qp_solver->orig_dims, c_ptr);
        c_ptr += ocp_qp_res_calculate_size(dims->qp_solver->orig_dims);

        // qp res ws
        work->nlp_work->qp_res_ws = ocp_qp_res_workspace_assign(dims->qp_solver->orig_dims, c_ptr);
        c_ptr += ocp_qp_res_workspace_calculate_size(dims->qp_solver->orig_dims);
    }

    assert((char *) work + ocp_nlp_sqp_wfqp_workspace_calculate_size(config, dims, opts) >= c_ptr);

    return;
}

void ocp_nlp_sqp_wfqp_work_get(void *config_, void *dims_, void *work_,
                          const char *field, void *return_value_)
{
    // ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_workspace *work = work_;

    if (!strcmp("nlp_work", field))
    {
        void **value = return_value_;
        *value = work->nlp_work;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_sqp_work_get\n", field);
        exit(1);
    }
}


/************************************************
 * termination criterion
 ************************************************/
static bool check_termination(int n_iter, ocp_nlp_dims *dims, ocp_nlp_res *nlp_res, ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_sqp_wfqp_opts *opts)
{
    // ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    // check for nans
    if (isnan(nlp_res->inf_norm_res_stat) || isnan(nlp_res->inf_norm_res_eq) ||
        isnan(nlp_res->inf_norm_res_ineq) || isnan(nlp_res->inf_norm_res_comp))
    {
        mem->nlp_mem->status = ACADOS_NAN_DETECTED;
        if (opts->nlp_opts->print_level > 0)
        {
            printf("Stopped: NaN detected in iterate.\n");
        }
        return true;
    }

    // check for maximum iterations
    if (!opts->eval_residual_at_max_iter && n_iter >= opts->max_iter)
    {
        mem->nlp_mem->status = ACADOS_MAXITER;
        if (opts->nlp_opts->print_level > 0)
        {
            printf("Stopped: Maximum Iterations Reached.\n");
        }
        return true;
    }

    // check if solved to tolerance
    if ((nlp_res->inf_norm_res_stat < opts->tol_stat) &&
        (nlp_res->inf_norm_res_eq < opts->tol_eq) &&
        (nlp_res->inf_norm_res_ineq < opts->tol_ineq) &&
        (nlp_res->inf_norm_res_comp < opts->tol_comp))
    {
        mem->nlp_mem->status = ACADOS_SUCCESS;
        if (opts->nlp_opts->print_level > 0)
        {
            printf("Optimal Solution found! Converged to KKT point.\n");
        }
        return true;
    }

    // check for small step
    if (opts->tol_min_step_norm > 0.0 && (n_iter > 0) && (mem->step_norm < opts->tol_min_step_norm))
    {
        if (opts->nlp_opts->print_level > 0)
        {
            if (nlp_res->inf_norm_res_eq < opts->tol_eq && nlp_res->inf_norm_res_ineq < opts->tol_ineq)
            {
                printf("Stopped: Converged to Feasible Point. Step size is < tol_eq.\n");
            }
            else
            {
                printf("Stopped: Converged to Infeasible Point. Step size is < tol_eq.\n");
            }
        }
        mem->nlp_mem->status = ACADOS_MINSTEP;
        return true;
    }

    // check for unbounded problem
    if (mem->nlp_mem->cost_value <= opts->tol_unbounded)
    {
        mem->nlp_mem->status = ACADOS_UNBOUNDED;
        if (opts->nlp_opts->print_level > 0)
        {
            printf("Stopped: Problem seems to be unbounded.\n");
        }
        return true;
    }

    // check for maximum iterations
    if (n_iter >= opts->max_iter)
    {
        mem->nlp_mem->status = ACADOS_MAXITER;
        if (opts->nlp_opts->print_level > 0)
        {
            printf("Stopped: Maximum Iterations Reached.\n");
        }
        return true;
    }

    return false;
}


/************************************************
 * functions
 ************************************************/
static void set_non_slacked_l1_penalties(ocp_nlp_config *config, ocp_nlp_dims *dims,
    ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_sqp_wfqp_memory *mem,
    ocp_nlp_workspace *work)
{
    int N = dims->N;
    // int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int *nns = mem->nns;
    ocp_qp_in *qp_in = mem->nlp_mem->qp_in;

    // be aware of rqz_QP = [r, q, zl_NLP, zl_QP, zu_NLP, zu_QP]
    for (int stage = 0; stage <= dims->N; stage++)
    {
        // zl_QP
        blasfeo_dvecse(nns[stage], mem->penalty_parameter, qp_in->rqz+stage, nu[stage]+nx[stage]+ns[stage]);
        // zu_QP
        blasfeo_dvecse(nns[stage], mem->penalty_parameter, qp_in->rqz+stage, nu[stage]+nx[stage]+2*ns[stage]+nns[stage]);
        // printf("qp_in->rqz %d\n", stage);
        // blasfeo_print_exp_tran_dvec(nu[stage] +nx[stage] + 2*(nns[stage]+ns[stage]), qp_in->rqz+stage, 0);
    }
}


static void set_non_slacked_l2_penalties(ocp_nlp_config *config, ocp_nlp_dims *dims,
    ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_sqp_wfqp_memory *mem,
    ocp_nlp_workspace *work)
{
    int N = dims->N;
    // int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int *nns = mem->nns;
    ocp_qp_in *qp_in = mem->nlp_mem->qp_in;

    // be aware of rqz_QP = [r, q, zl_NLP, zl_QP, zu_NLP, zu_QP]
    for (int stage = 0; stage <= dims->N; stage++)
    {
        // zu_NLP shift back
        blasfeo_dveccp(ns[stage], qp_in->Z+stage, ns[stage], qp_in->Z+stage, ns[stage]+nns[stage]);
        // zl_QP
        blasfeo_dvecse(nns[stage], 0.0, qp_in->Z+stage, ns[stage]);
        // zu_QP
        blasfeo_dvecse(nns[stage], 0.0, qp_in->Z+stage, 2*ns[stage]+nns[stage]);

        // printf("qp_in->Z %d\n", stage);
        // blasfeo_print_exp_tran_dvec(2*(nns[stage]+ns[stage]), qp_in->Z+stage, 0);
    }
    // print_ocp_qp_in(qp_in);
}



// TODO: work this out!
// signature: need solver specific memory
/*
calculates new iterate or trial iterate in 'out_destination' with step 'mem->qp_out',
step size 'alpha', and current iterate 'out_start'.
 */
void ocp_nlp_update_variables_sqp_wfqp(void *config_, void *dims_,
            void *in_, void *out_, void *opts_, void *mem_,
            void *work_, void *out_destination_,
            void *solver_mem, double alpha, bool full_step_dual)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_out *out_start = out_;
    ocp_nlp_memory *nlp_mem = mem_;
    ocp_nlp_sqp_wfqp_memory *mem = solver_mem;
    
    ocp_nlp_out *out_destination = out_destination_;
    // solver_mem is not used in this function, but needed for DDP
    // the function is used in the config->globalization->step_update
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;
    int *nz = dims->nz;
    int *ns = dims->ns;
    int *nns = mem->nns;
    int n_nominal_ineq_nlp;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // step in primal variables
        blasfeo_daxpy(nx[i]+nu[i]+ns[i], alpha, nlp_mem->qp_out->ux + i, 0, out_start->ux + i, 0, out_destination->ux + i, 0);
        blasfeo_daxpy(ns[i], alpha, nlp_mem->qp_out->ux + i, nx[i]+nu[i]+ns[i]+nns[i], out_start->ux + i, nx[i]+nu[i]+ns[i], out_destination->ux + i, nx[i]+nu[i]+ns[i]);

        // update dual variables
        n_nominal_ineq_nlp = 2*dims->ni[i] -2*ns[i]; //additional slacks are not counted in dims->ni[i]
        // Assuming a constraint order
        // [lbu, ubu, lbx, ubx, lbg, ubg, lbh, ubh, lbs, ubs]
        // n_nominal_ineq_nlp = [lbu, ubu, lbx, ubx, lbg, ubg, lbh, ubh]
        if (full_step_dual)
        {
            // printf("at stage: %d, ni[i]= %d, ns[i]=%d, nns[i]= %d\n", i, dims->ni[i], ns[i], mem->nns[i]);
            // TODO: split ni_qp = n_ineq_nominal_nlp + ns_nlp + nns;
            blasfeo_dveccp(n_nominal_ineq_nlp+ns[i], nlp_mem->qp_out->lam+i, 0, out_destination->lam+i, 0);
            blasfeo_dveccp(ns[i], nlp_mem->qp_out->lam+i, n_nominal_ineq_nlp+ns[i]+mem->nns[i], out_destination->lam+i, n_nominal_ineq_nlp+ns[i]+mem->nns[i]);
            if (i < N)
            {
                blasfeo_dveccp(nx[i+1], nlp_mem->qp_out->pi+i, 0, out_destination->pi+i, 0);
            }
        }
        else
        {
            // TODO: split ni_qp = n_ineq_nominal_nlp + ns_nlp + nns;
            // update duals with alpha step
            blasfeo_daxpby(n_nominal_ineq_nlp+ns[i], 1.0-alpha, out_start->lam+i, 0, alpha, nlp_mem->qp_out->lam+i, 0, out_destination->lam+i, 0);
            blasfeo_daxpby(ns[i], 1.0-alpha, out_start->lam+i, n_nominal_ineq_nlp+ns[i]+mem->nns[i], alpha, nlp_mem->qp_out->lam+i, n_nominal_ineq_nlp+ns[i]+mem->nns[i], out_destination->lam+i, n_nominal_ineq_nlp+ns[i]+mem->nns[i]);
            
            // blasfeo_dvecsc(2*ni[i], 1.0-alpha, out->lam+i, 0);
            // blasfeo_daxpy(2*ni[i], alpha, nlp_mem->qp_out->lam+i, 0, out->lam+i, 0, out->lam+i, 0);
            if (i < N)
            {
                // blasfeo_dvecsc(nx[i+1], 1.0-alpha, out->pi+i, 0);
                // blasfeo_daxpy(nx[i+1], alpha, nlp_mem->qp_out->pi+i, 0, out->pi+i, 0, out->pi+i, 0);
                blasfeo_daxpby(nx[i+1], 1.0-alpha, out_start->pi+i, 0, alpha, nlp_mem->qp_out->pi+i, 0, out_destination->pi+i, 0);
            }
        }

        // linear update of algebraic variables using state and input sensitivity
        if (i < N)
        {
            // out->z = nlp_mem->z_alg + alpha * dzdux * qp_out->ux
            blasfeo_dgemv_t(nu[i]+nx[i], nz[i], alpha, nlp_mem->dzduxt+i, 0, 0,
                    nlp_mem->qp_out->ux+i, 0, 1.0, nlp_mem->z_alg+i, 0, out_destination->z+i, 0);
        }
    }
}


static void print_indices(
    ocp_nlp_dims *dims,
    ocp_nlp_sqp_wfqp_memory *mem
)
{
    // DEBUG:
    for (int stage = 0; stage <= dims->N; stage++)
    {
        int nns = mem->nns[stage];
        printf("got idxns at stage %d\n", stage);
        for (int i=0; i<nns; i++)
            printf("%d ", mem->idxns[stage][i]);
        printf("\n");
    }
}


void ocp_nlp_approximate_qp_vectors_sqp_wfqp(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
    ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_workspace *work)
{
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int *ni = dims->ni;
    int *nns = mem->nns;


#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // g
        // standard
        // blasfeo_dveccp(nv[i], nlp_mem->cost_grad + i, 0, nlp_mem->qp_in->rqz + i, 0);
        blasfeo_dveccp(nx[i]+nu[i]+ns[i], nlp_mem->cost_grad + i, 0, nlp_mem->qp_in->rqz + i, 0);
        blasfeo_dveccp(ns[i], nlp_mem->cost_grad + i, nx[i]+nu[i]+ns[i], nlp_mem->qp_in->rqz + i, nx[i]+nu[i]+ns[i]+nns[i]);

        // b
        if (i < N)
            blasfeo_dveccp(nx[i + 1], nlp_mem->dyn_fun + i, 0, nlp_mem->qp_in->b + i, 0);

        // evaluate constraint residuals
        config->constraints[i]->update_qp_vectors(config->constraints[i], dims->constraints[i],
            in->constraints[i], opts->constraints[i], nlp_mem->constraints[i], work->constraints[i]);

        // TODO: ni!
        // copy ineq function value into nlp mem, then into QP
        struct blasfeo_dvec *ineq_fun = config->constraints[i]->memory_get_fun_ptr(nlp_mem->constraints[i]);
        blasfeo_dveccp(2 * ni[i], ineq_fun, 0, nlp_mem->ineq_fun + i, 0);

        // d
        blasfeo_dveccp(2 * ni[i], nlp_mem->ineq_fun + i, 0, nlp_mem->qp_in->d + i, 0);
    }
}


// MAIN OPTIMIZATION ROUTINE
int ocp_nlp_sqp_wfqp(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_)
{
    acados_timer timer0, timer1;
    acados_tic(&timer0);

    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_sqp_wfqp_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_res *nlp_res = nlp_mem->nlp_res;
    ocp_nlp_timings *nlp_timings = nlp_mem->nlp_timings;

    ocp_nlp_sqp_wfqp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    ocp_qp_in *qp_in = nlp_mem->qp_in;
    ocp_qp_out *qp_out = nlp_mem->qp_out;

    // zero timers
    ocp_nlp_timings_reset(nlp_timings);

    int qp_status = 0;
    int qp_iter = 0;
    mem->alpha = 0.0;
    mem->step_norm = 0.0;
    mem->nlp_mem->status = ACADOS_SUCCESS;

    // TODO(@david):
    mem->penalty_parameter = 42*1e5;

#if defined(ACADOS_WITH_OPENMP)
    // backup number of threads
    int num_threads_bkp = omp_get_num_threads();
    // set number of threads
    omp_set_num_threads(opts->nlp_opts->num_threads);
#endif

    ocp_nlp_initialize_submodules(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
    set_non_slacked_l2_penalties(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work);

    /************************************************
     * main sqp loop
     ************************************************/
    int sqp_iter = 0;
    double prev_levenberg_marquardt = 0.0;
    for (; sqp_iter <= opts->max_iter; sqp_iter++) // <= needed such that after last iteration KKT residuals are checked before max_iter is thrown.
    {
        // We always evaluate the residuals until the last iteration
        // If the option "eval_residual_at_max_iter" is set, we also
        // evaluate the residuals after the last iteration.
        if (sqp_iter != opts->max_iter || opts->eval_residual_at_max_iter)
        {
            /* Prepare the QP data */
            // linearize NLP and update QP matrices
            acados_tic(&timer1);
            ocp_nlp_approximate_qp_matrices(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);

            if (nlp_opts->with_adaptive_levenberg_marquardt || config->globalization->needs_objective_value() == 1)
            {
                ocp_nlp_get_cost_value_from_submodules(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
            }
            ocp_nlp_add_levenberg_marquardt_term(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, mem->alpha, sqp_iter);

            // update QP rhs for SQP (step prim var, abs dual var)
            ocp_nlp_approximate_qp_vectors_sqp_wfqp(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work);
            nlp_timings->time_lin += acados_toc(&timer1);

            // Set the penalties in slacked problem
            set_non_slacked_l1_penalties(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work);
            // compute nlp residuals
            ocp_nlp_res_compute(dims, nlp_in, nlp_out, nlp_res, nlp_mem);
            ocp_nlp_res_get_inf_norm(nlp_res, &nlp_out->inf_norm_res);
        }

        // Initialize the memory for different globalization strategies
        if (sqp_iter == 0)
        {
            config->globalization->initialize_memory(config, dims, nlp_mem, nlp_opts);
        }

        // save statistics
        if (sqp_iter < mem->stat_m)
        {
            mem->stat[mem->stat_n*sqp_iter+0] = nlp_res->inf_norm_res_stat;
            mem->stat[mem->stat_n*sqp_iter+1] = nlp_res->inf_norm_res_eq;
            mem->stat[mem->stat_n*sqp_iter+2] = nlp_res->inf_norm_res_ineq;
            mem->stat[mem->stat_n*sqp_iter+3] = nlp_res->inf_norm_res_comp;
        }

        // Output
        if (nlp_opts->print_level > 0)
        {
            config->globalization->print_iteration(nlp_mem->cost_value,
                                                   sqp_iter,
                                                   nlp_res,
                                                   mem->step_norm,
                                                   prev_levenberg_marquardt,
                                                   qp_status,
                                                   qp_iter,
                                                   nlp_opts,
                                                   nlp_mem->globalization);
        }
        prev_levenberg_marquardt = nlp_opts->levenberg_marquardt;

        // regularize Hessian
        // NOTE: this is done before termination, such that we can get the QP at the stationary point that is actually solved, if we exit with success.
        acados_tic(&timer1);
        config->regularize->regularize(config->regularize, dims->regularize,
                                               nlp_opts->regularize, nlp_mem->regularize_mem);
        nlp_timings->time_reg += acados_toc(&timer1);

        // Termination
        if (check_termination(sqp_iter, dims, nlp_res, mem, opts))
        {
#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif
            nlp_mem->iter = sqp_iter;
            nlp_timings->time_tot = acados_toc(&timer0);
            return mem->nlp_mem->status;
        }

        /* solve QP */
        // (typically) no warm start at first iteration
        if (sqp_iter == 0 && !opts->warm_start_first_qp)
        {
            int tmp_int = 0;
            qp_solver->opts_set(qp_solver, nlp_opts->qp_solver_opts, "warm_start", &tmp_int);
        }
        // Show input to QP
        if (nlp_opts->print_level > 3)
        {
            printf("\n\nSQP: ocp_qp_in at iteration %d\n", sqp_iter);
            print_ocp_qp_in(qp_in);
        }

#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
        ocp_nlp_dump_qp_in_to_file(qp_in, sqp_iter, 0);
#endif

        // printf("QP in before solve\n");
        // print_ocp_qp_in(qp_in);
        // print_indices(dims, mem);
        // exit(1);

        qp_status = ocp_nlp_solve_qp_and_correct_dual(config, dims, nlp_opts, nlp_mem, nlp_work, false);

        // restore default warm start
        if (sqp_iter==0)
        {
            qp_solver->opts_set(qp_solver, nlp_opts->qp_solver_opts, "warm_start", &opts->qp_warm_start);
        }

        if (nlp_opts->print_level > 3)
        {
            printf("\n\nSQP: ocp_qp_out at iteration %d\n", sqp_iter);
            print_ocp_qp_out(qp_out);
        }

#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
        ocp_nlp_dump_qp_out_to_file(qp_out, sqp_iter, 0);
#endif

        qp_info *qp_info_;
        ocp_qp_out_get(qp_out, "qp_info", &qp_info_);
        qp_iter = qp_info_->num_iter;

        // save statistics of last qp solver call
        if (sqp_iter+1 < mem->stat_m)
        {
            mem->stat[mem->stat_n*(sqp_iter+1)+4] = qp_status;
            mem->stat[mem->stat_n*(sqp_iter+1)+5] = qp_iter;
        }

        // compute external QP residuals (for debugging)
        if (nlp_opts->ext_qp_res)
        {
            ocp_qp_res_compute(qp_in, qp_out, nlp_work->qp_res, nlp_work->qp_res_ws);
            if (sqp_iter+1 < mem->stat_m)
                ocp_qp_res_compute_nrm_inf(nlp_work->qp_res, mem->stat+(mem->stat_n*(sqp_iter+1)+7));
        }

        // exit conditions on QP status
        if ((qp_status!=ACADOS_SUCCESS) & (qp_status!=ACADOS_MAXITER))
        {
            if (nlp_opts->print_level > 0)
            {
                printf("%i\t%e\t%e\t%e\t%e.\n", sqp_iter, nlp_res->inf_norm_res_stat,
                    nlp_res->inf_norm_res_eq, nlp_res->inf_norm_res_ineq,
                    nlp_res->inf_norm_res_comp );
                printf("\n\n");
            }
            // increment sqp_iter to return full statistics and improve output below.
            sqp_iter++;

#ifndef ACADOS_SILENT
            printf("\nQP solver returned error status %d in SQP iteration %d, QP iteration %d.\n",
                   qp_status, sqp_iter, qp_iter);
#endif
#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif
            if (nlp_opts->print_level > 1)
            {
                printf("\n Failed to solve the following QP:\n");
                if (nlp_opts->print_level)
                    print_ocp_qp_in(qp_in);
            }

            mem->nlp_mem->status = ACADOS_QP_FAILURE;
            nlp_mem->iter = sqp_iter;
            nlp_timings->time_tot = acados_toc(&timer0);

            return mem->nlp_mem->status;
        }

        // Calculate optimal QP objective (needed for globalization)
        if (config->globalization->needs_qp_objective_value() == 1)
        {
            nlp_mem->qp_cost_value = ocp_nlp_compute_qp_objective_value(dims, qp_in, qp_out, nlp_work);
        }

        // Compute the step norm
        if (opts->tol_min_step_norm > 0.0 || nlp_opts->log_primal_step_norm)
        {
            mem->step_norm = ocp_qp_out_compute_primal_nrm_inf(nlp_mem->qp_out);
            if (nlp_opts->log_primal_step_norm)
                mem->primal_step_norm[sqp_iter] = mem->step_norm;
        }
        /* end solve QP */

        /* globalization */
        // NOTE on timings: currently all within globalization is accounted for within time_glob.
        //   QP solver times could be also attributed there alternatively. Cleanest would be to save them seperately.
        acados_tic(&timer1);

        int globalization_status;
        globalization_status = config->globalization->find_acceptable_iterate(config, dims, nlp_in, nlp_out, nlp_mem, mem, nlp_work, nlp_opts, &mem->alpha);

        if (globalization_status != ACADOS_SUCCESS)
        {
            if (nlp_opts->print_level > 1)
            {
                printf("\nFailure in globalization, got status %d!\n", globalization_status);
            }
            mem->nlp_mem->status = globalization_status;
            nlp_mem->iter = sqp_iter;
            nlp_timings->time_tot = acados_toc(&timer0);
#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif
            return mem->nlp_mem->status;
        }
        mem->stat[mem->stat_n*(sqp_iter+1)+6] = mem->alpha;
        nlp_timings->time_glob += acados_toc(&timer1);

    }  // end SQP loop

    if (nlp_opts->print_level > 0)
    {
        printf("Warning: The solver should never reach this part of the function!\n");
    }
#if defined(ACADOS_WITH_OPENMP)
    // restore number of threads
    omp_set_num_threads(num_threads_bkp);
#endif
    return mem->nlp_mem->status;
}

double ocp_nlp_sqp_wfqp_compute_qp_objective_value(ocp_nlp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_nlp_workspace *nlp_work)
{
    // Compute the QP objective function value
    double qp_cost = 0.0;
    int i, nux, ns;
    int N = dims->N;
    // Sum over stages 0 to N
    for (i = 0; i <= N; i++)
    {
        nux = dims->nx[i] + dims->nu[i];
        ns = dims->ns[i];
        // Calculate 0.5 * d.T H d
        blasfeo_dsymv_l(nux, 0.5, &qp_in->RSQrq[i], 0, 0, &qp_out->ux[i], 0,
                        0.0, &qp_out->ux[i], 0, &nlp_work->tmp_nv, 0);
        qp_cost += blasfeo_ddot(nux, &qp_out->ux[i], 0, &nlp_work->tmp_nv, 0);

        // slack QP objective value, compare to computation in cost modules;
        // tmp_nv = 2 * z + Z .* slack;
        blasfeo_dveccpsc(2*ns, 2.0, &qp_out->ux[i], nux, &nlp_work->tmp_nv, 0);
        blasfeo_dvecmulacc(2*ns, &qp_in->Z[i], 0, &qp_out->ux[i], nux, &nlp_work->tmp_nv, 0);
        // qp_cost += .5 * (tmp_nv .* slack)
        qp_cost += 0.5 * blasfeo_ddot(2*ns, &nlp_work->tmp_nv, 0, &qp_out->ux[i], nux);
        // Calculate g.T d
        qp_cost += blasfeo_ddot(nux, &qp_out->ux[i], 0, &qp_in->rqz[i], 0);

        // Calculate gradient of slacks
        qp_cost += blasfeo_ddot(2 * ns, &qp_out->ux[i], nux, &qp_in->rqz[i], nux);
    }
    return qp_cost;
}

void ocp_nlp_sqp_wfqp_memory_reset_qp_solver(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
    void *opts_, void *mem_, void *work_)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_opts *opts = opts_;
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_sqp_wfqp_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_sqp_wfqp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    // printf("in ocp_nlp_sqp_wfqp_memory_reset_qp_solver\n\n");
    config->qp_solver->memory_reset(qp_solver, dims->qp_solver,
        nlp_mem->qp_in, nlp_mem->qp_out, opts->nlp_opts->qp_solver_opts,
        nlp_mem->qp_solver_mem, nlp_work->qp_work);
}


int ocp_nlp_sqp_wfqp_precompute(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_opts *opts = opts_;
    ocp_nlp_sqp_wfqp_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    nlp_mem->workspace_size = ocp_nlp_workspace_calculate_size(config, dims, opts->nlp_opts);

    ocp_nlp_sqp_wfqp_workspace *work = work_;
    ocp_nlp_sqp_wfqp_cast_workspace(config, dims, opts, mem, work);
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    // create indices
    int ni, ns, nns, nbu, nbx;
    int *idxs = nlp_work->tmp_nins;
    for (int stage = 0; stage <= dims->N; stage++)
    {
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage], "nbu", &nbu);
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage], "nbx", &nbx);
        ns = dims->ns[stage];
        ni = dims->ni[stage];
        nns = mem->nns[stage];
        config->constraints[stage]->model_get(config->constraints[stage], dims->constraints[stage], nlp_in->constraints[stage], "idxs", idxs);

        // DEBUG:
        printf("stage %d: ni %d ns %d nns %d\n", stage, ni, ns, nns);
        printf("got idxs at stage %d\n", stage);
        for (int i=0; i<ns; i++)
            printf("%d ", idxs[i]);
        printf("\n");

        int ins = 0;
        bool i_slacked;
        // create idxns
        for (int i=0; i<ni-ns; i++)
        {
            i_slacked = false;
            for (int j=0; j<ns; j++)
            {
                if (idxs[j] == i)
                {
                    i_slacked = true;
                    break;
                }
            }
            // slack constraints that are not slacked yet and not u-bounds
            if (!i_slacked && ((stage>0 && i >= nbu) || (stage==0 && i >= nbu+nbx)))
            {
                mem->idxns[stage][ins] = i;
                ins++;
            }
        }
        // DEBUG:
        printf("got idxns at stage %d\n", stage);
        for (int i=0; i<nns; i++)
            printf("%d ", mem->idxns[stage][i]);
        if (ins != nns)
        {
            printf("ins %d != nns %d", ins, nns);
        }
        printf("\n");
    }


    // set idxs_rev of QP:
    for (int stage = 0; stage <= dims->N; stage++)
    {
        config->constraints[stage]->model_get(config->constraints[stage], dims->constraints[stage], nlp_in->constraints[stage], "idxs", idxs);
        int *idxs_rev = nlp_mem->qp_in->idxs_rev[stage];
        int *idxns = mem->idxns[stage];
        for (int i=0; i<nns; i++)
        {
            idxs_rev[idxns[i]] = i+dims->ns[stage];
        }
    }

    return ocp_nlp_precompute_common(config, dims, nlp_in, nlp_out, opts->nlp_opts, nlp_mem, nlp_work);
}



void ocp_nlp_sqp_wfqp_eval_param_sens(void *config_, void *dims_, void *opts_, void *mem_, void *work_,
                                 char *field, int stage, int index, void *sens_nlp_out_)
{
    acados_timer timer0;
    acados_tic(&timer0);

    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_opts *opts = opts_;
    ocp_nlp_sqp_wfqp_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_nlp_out *sens_nlp_out = sens_nlp_out_;

    ocp_nlp_sqp_wfqp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    ocp_nlp_common_eval_param_sens(config, dims, opts->nlp_opts, nlp_mem, nlp_work,
                                 field, stage, index, sens_nlp_out);

    nlp_mem->nlp_timings->time_solution_sensitivities = acados_toc(&timer0);

    return;
}


void ocp_nlp_sqp_wfqp_eval_lagr_grad_p(void *config_, void *dims_, void *nlp_in_, void *opts_, void *mem_, void *work_,
                                 const char *field, void *grad_p)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_opts *opts = opts_;
    ocp_nlp_sqp_wfqp_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    ocp_nlp_in *nlp_in = nlp_in_;

    ocp_nlp_sqp_wfqp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    ocp_nlp_common_eval_lagr_grad_p(config, dims, nlp_in, opts->nlp_opts, nlp_mem, nlp_work,
                                 field, grad_p);

    return;
}


void ocp_nlp_sqp_wfqp_get(void *config_, void *dims_, void *mem_, const char *field, void *return_value_)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_sqp_wfqp_memory *mem = mem_;


    char module[MAX_STR_LEN];
    char *ptr_module = NULL;
    int module_length = 0;

    // extract module name
    char *char_ = strchr(field, '_');
    if (char_!=NULL)
    {
        module_length = char_-field;
        for (int ii=0; ii<module_length; ii++)
            module[ii] = field[ii];
        module[module_length] = '\0'; // add end of string
        ptr_module = module;
    }

    if ( ptr_module!=NULL && (!strcmp(ptr_module, "time")) )
    {
        // call timings getter
        ocp_nlp_timings_get(config, mem->nlp_mem->nlp_timings, field, return_value_);
    }
    else if (!strcmp("stat", field))
    {
        double **value = return_value_;
        *value = mem->stat;
    }
    else if (!strcmp("primal_step_norm", field))
    {
        if (mem->primal_step_norm == NULL)
        {
            printf("\nerror: options log_primal_step_norm was not set\n");
            exit(1);
        }
        else
        {
            double *value = return_value_;
            for (int ii=0; ii<mem->nlp_mem->iter; ii++)
            {
                value[ii] = mem->primal_step_norm[ii];
            }
        }
    }
    else if (!strcmp("statistics", field))
    {
        int n_row = mem->stat_m<mem->nlp_mem->iter+1 ? mem->stat_m : mem->nlp_mem->iter+1;
        double *value = return_value_;
        for (int ii=0; ii<n_row; ii++)
        {
            value[ii+0] = ii;
            for (int jj=0; jj<mem->stat_n; jj++)
                value[ii+(jj+1)*n_row] = mem->stat[jj+ii*mem->stat_n];
        }
    }
    else if (!strcmp("stat_m", field))
    {
        int *value = return_value_;
        *value = mem->stat_m;
    }
    else if (!strcmp("stat_n", field))
    {
        int *value = return_value_;
        *value = mem->stat_n;
    }
    else if (!strcmp("qp_xcond_dims", field))
    {
        void **value = return_value_;
        *value = dims->qp_solver->xcond_dims;
    }
    else
    {
        ocp_nlp_memory_get(config, mem->nlp_mem, field, return_value_);
    }
}


void ocp_nlp_sqp_wfqp_terminate(void *config_, void *mem_, void *work_)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_memory *mem = mem_;
    ocp_nlp_sqp_wfqp_workspace *work = work_;

    config->qp_solver->terminate(config->qp_solver, mem->nlp_mem->qp_solver_mem, work->nlp_work->qp_work);
}

bool ocp_nlp_sqp_wfqp_is_real_time_algorithm()
{
    return false;
}

void ocp_nlp_sqp_wfqp_config_initialize_default(void *config_)
{
    ocp_nlp_config *config = (ocp_nlp_config *) config_;

    config->with_feasible_qp = 1;

    config->opts_calculate_size = &ocp_nlp_sqp_wfqp_opts_calculate_size;
    config->opts_assign = &ocp_nlp_sqp_wfqp_opts_assign;
    config->opts_initialize_default = &ocp_nlp_sqp_wfqp_opts_initialize_default;
    config->opts_update = &ocp_nlp_sqp_wfqp_opts_update;
    config->opts_set = &ocp_nlp_sqp_wfqp_opts_set;
    config->opts_set_at_stage = &ocp_nlp_sqp_wfqp_opts_set_at_stage;
    config->memory_calculate_size = &ocp_nlp_sqp_wfqp_memory_calculate_size;
    config->memory_assign = &ocp_nlp_sqp_wfqp_memory_assign;
    config->workspace_calculate_size = &ocp_nlp_sqp_wfqp_workspace_calculate_size;
    config->evaluate = &ocp_nlp_sqp_wfqp;
    config->memory_reset_qp_solver = &ocp_nlp_sqp_wfqp_memory_reset_qp_solver;
    config->eval_param_sens = &ocp_nlp_sqp_wfqp_eval_param_sens;
    config->eval_lagr_grad_p = &ocp_nlp_sqp_wfqp_eval_lagr_grad_p;
    config->config_initialize_default = &ocp_nlp_sqp_wfqp_config_initialize_default;
    config->precompute = &ocp_nlp_sqp_wfqp_precompute;
    config->get = &ocp_nlp_sqp_wfqp_get;
    config->opts_get = &ocp_nlp_sqp_wfqp_opts_get;
    config->work_get = &ocp_nlp_sqp_wfqp_work_get;
    config->terminate = &ocp_nlp_sqp_wfqp_terminate;
    config->step_update = &ocp_nlp_update_variables_sqp_wfqp;
    config->is_real_time_algorithm = &ocp_nlp_sqp_wfqp_is_real_time_algorithm;

    return;
}