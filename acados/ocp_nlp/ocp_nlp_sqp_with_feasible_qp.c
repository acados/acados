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
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
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
    opts->initial_objective_multiplier = 1e0;
    opts->sufficient_l1_inf_reduction = 1e-1;

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
        else if (!strcmp(field, "initial_objective_multiplier"))
        {
            double* initial_objective_multiplier = (double *) value;
            opts->initial_objective_multiplier = *initial_objective_multiplier;
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
    int N = dims->N;

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
    size += (N + 1) * sizeof(int *);
    int nns, nsbu, nbu, nsbx, nbx, n_nominal_ineq_nlp;
    for (int stage = 0; stage <= N; stage++)
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
        // idxns
        size += nns * sizeof(int);
        // slacks_not_in_original_nlp
        size += blasfeo_memsize_dvec(2*nns);
        // Z_cost_module
        size += blasfeo_memsize_dvec(2*dims->ns[stage]);

        // RSQ_cost, RSQ_constr
        size += 2*blasfeo_memsize_dmat(dims->nx[stage]+dims->nu[stage], dims->nx[stage]+dims->nu[stage]);
    }
    // nns
    size += (N+1) * sizeof(int);

    // slacks_not_in_original_nlp
    size += (N + 1) * sizeof(struct blasfeo_dvec);
    // Z_cost_module
    size += (N + 1) * sizeof(struct blasfeo_dvec);
    // RSQ_cost, RSQ_constr
    size += 2*(N + 1) * sizeof(struct blasfeo_dmat);

    size += 3*8;  // align
    size += 64;  // blasfeo_mem align

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

    // slacks_not_in_original_nlp
    assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->slacks_not_in_original_nlp, &c_ptr);
    // Z_cost_module
    assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->Z_cost_module, &c_ptr);

    // RSQ_cost
    assign_and_advance_blasfeo_dmat_structs(N + 1, &mem->RSQ_cost, &c_ptr);
    // RSQ_constr
    assign_and_advance_blasfeo_dmat_structs(N + 1, &mem->RSQ_constr, &c_ptr);

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

    int nsbu, nbu, nsbx, nbx, n_nominal_ineq_nlp;
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

    // blasfeo_mem align
    align_char_to(64, &c_ptr);
    // matrices first: RSQ
    for (int i = 0; i <= N; ++i)
    {
        assign_and_advance_blasfeo_dmat_mem(dims->nx[i]+dims->nu[i], dims->nx[i]+dims->nu[i], mem->RSQ_cost + i, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(dims->nx[i]+dims->nu[i], dims->nx[i]+dims->nu[i], mem->RSQ_constr + i, &c_ptr);
    }
    // blasfeo_dvec
    // slacks_not_in_original_nlp
    for (int i = 0; i <= N; ++i)
    {
        assign_and_advance_blasfeo_dvec_mem(2*mem->nns[i], mem->slacks_not_in_original_nlp + i, &c_ptr);
    }
    // Z_cost_module
    for (int i = 0; i <= N; ++i)
    {
        assign_and_advance_blasfeo_dvec_mem(2*dims->ns[i], mem->Z_cost_module + i, &c_ptr);
    }
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

static double get_multiplier_norm_inf(ocp_nlp_dims *dims, ocp_qp_out *qp_out)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *ns = dims->ns;
    double tmp0 = 0.0;

    for (int i = 0; i < N; i++)
    {
        for (int j=0; j<nx[i+1]; j++)
        {
            // abs(lambda) (LW)
            tmp0 = fmax(tmp0, fabs(BLASFEO_DVECEL(qp_out->pi+i, j)));
        }
    }
    for (int i = 0; i <= N; i++)
    {
        int n_nominal_ineq_nlp = dims->ni[i] - ns[i];
        int two_n_nominal_ineq_nlp = 2*n_nominal_ineq_nlp;
        for (int j=0; j<two_n_nominal_ineq_nlp+ns[i]; j++)
        {
            // mu (LW)
            tmp0 = fmax(tmp0, BLASFEO_DVECEL(qp_out->lam+i, j));
        }
        for (int j=0; j<ns[i]; j++)
        {
            // mu (LW)
            tmp0 = fmax(tmp0, BLASFEO_DVECEL(qp_out->lam+i, two_n_nominal_ineq_nlp+ns[i]+j));
        }
    }
    return tmp0;
}

static void scale_multiplier(ocp_nlp_dims *dims, ocp_nlp_memory* nlp_mem, ocp_qp_out *qp_out)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *ns = dims->ns;

    for (int i = 0; i < N; i++)
    {
        blasfeo_dvecsc(nx[i+1], 1.0/nlp_mem->objective_multiplier, &qp_out->pi[i], 0);
    }
    for (int i = 0; i <= N; i++)
    {
        int n_nominal_ineq_nlp = dims->ni[i] - ns[i];
        int two_n_nominal_ineq_nlp = 2*n_nominal_ineq_nlp;
        blasfeo_dvecsc(two_n_nominal_ineq_nlp+ns[i], 1.0/nlp_mem->objective_multiplier, qp_out->lam+i, 0);
        blasfeo_dvecsc(ns[i], 1.0/nlp_mem->objective_multiplier, qp_out->lam+i, two_n_nominal_ineq_nlp+ns[i]);
    }
}

static double calculate_predicted_l1_inf_reduction(ocp_nlp_sqp_wfqp_opts* opts, double current_infeasibility, double qp_infeasibility)
{
    if (current_infeasibility < fmin(opts->tol_ineq, opts->tol_eq))
    {
        return 0.0;
    }
    else
    {
        return current_infeasibility - qp_infeasibility;
    }
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
/*
This function calculates the convex interpolation factor
*/
static double calculate_search_direction_interpolation_factor(ocp_nlp_sqp_wfqp_opts* opts,
                                                              double pred_l1_inf_QP_feasibility,
                                                              double pred_l1_inf_QP_optimality,
                                                              double l1_inf_QP_optimality,
                                                              double l1_inf_QP_feasibility)
{
    double kappa;

    if (pred_l1_inf_QP_optimality >= opts->sufficient_l1_inf_reduction * pred_l1_inf_QP_feasibility)
    {
        kappa = 1.0;
    }    
    else
    {
        kappa = fmin(1.0, ((1-opts->sufficient_l1_inf_reduction)*pred_l1_inf_QP_feasibility)/(l1_inf_QP_optimality - l1_inf_QP_feasibility));
    }
    // We have a convex combination, therefore kappa in [0,1]
    // assert(kappa >= 0.0 && kappa <= 1.0);
    if (kappa < 0.0)
    {
        kappa = 0.0;
    }
    else if (kappa > 1.0)
    {
        kappa = 1.0;
    }
    printf("Kappa: %.4e\n", kappa);
    return kappa;
}

static void setup_search_direction(ocp_nlp_sqp_wfqp_memory* mem, ocp_nlp_dims* dims, ocp_qp_out* prediction_direction,
                                ocp_qp_out* steering_direction, ocp_qp_out* search_direction, double kappa)
{
    // steering direction is solution of feasibility problem
    // prediction direction is solution of optimality problem
    for (int i = 0; i <= dims->N; i++)
    {
        int dim = dims->nx[i]+dims->nu[i]+2*dims->ns[i]+2*mem->nns[i];
        blasfeo_daxpby(dim, kappa, &prediction_direction->ux[i], 0, 1-kappa, &steering_direction->ux[i], 0, &search_direction->ux[i], 0);
    }
}

/*
Given the current iterate and the value of the constraints at the current
iterate, this function calculates the minimal l1-slack variables which are not
included in the original NLP
*/
static void set_slack_variable_values(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_sqp_wfqp_memory *mem)
{
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    int stage, nns, i;
    double tmp;

    for (stage = 0; stage <= dims->N; stage++)
    {
        nns = mem->nns[stage];
        for (i=0; i<nns; i++)
        {
            tmp = BLASFEO_DVECEL(nlp_mem->ineq_fun+stage, mem->idxns[stage][i]);
            blasfeo_dvecin1(fmax(0.0, tmp), &mem->slacks_not_in_original_nlp[stage], i);
        }
    }
}

/*
This function calculates the l1 infeasibility based on the slack variable values
*/
static double get_slacked_qp_l1_infeasibility(ocp_nlp_dims *dims, ocp_nlp_sqp_wfqp_memory *mem, ocp_qp_out *qp_out)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int *nns = mem->nns;

    double l1_inf = 0.0;
    int i, j;
    double tmp1, tmp2;

    for (i = 0; i <= N; i++)
    {
        for (j=0; j<nns[i]; j++)
        {
            // Add lb slack
            tmp1 = BLASFEO_DVECEL(qp_out->ux + i, nx[i]+nu[i]+ns[i] + j);
            l1_inf += fmax(0.0, tmp1);
            // Add ub slack
            tmp2 = BLASFEO_DVECEL(qp_out->ux + i, nx[i]+nu[i]+2*ns[i]+nns[i] + j);
            l1_inf += fmax(0.0, tmp2);
        }
    }

    return l1_inf;
}

/*
This function calculates the l1 infeasibility by calculating the matrix vector product of the 
constraints
*/
static double manually_calculate_slacked_qp_l1_infeasibility(ocp_nlp_dims *dims, ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_sqp_wfqp_workspace *work, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nns = mem->nns;

    int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng; //number of general two sided constraints

    double l1_inf = 0.0;
    int i, j;
    double tmp, tmp_bound, mask_value;

    for (i = 0; i <= N; i++)
    {
        // bounds on states and controls
        for (j=0; j<nns[i]; ++j)
        {
            int index = mem->idxns[i][j];

            // get entry in the middle of inequality equation:
            // simple bounds
            if (index < nb[i])
            {
                tmp = BLASFEO_DVECEL(qp_out->ux+i, index);
            }
            // linear constraints
            else
            {
                // general linear / linearized!
                // tmp_ni = D * u + C * x
                // Calculate the product
                blasfeo_dgemv_t(nu[i]+nx[i], 1, 1.0, qp_in->DCt+i, index-nb[i], 0, qp_out->ux+i, 0,
                        0.0, qp_in->d+i, nb[i], &work->nlp_work->tmp_ni+i, 0);
                tmp = BLASFEO_DVECEL(&work->nlp_work->tmp_ni+i, 0);
            }
            
            // check lower bounds
            mask_value = BLASFEO_DVECEL(qp_in->d_mask+i, index);
            if (mask_value == 1.0)
            {
                tmp_bound = BLASFEO_DVECEL(qp_in->d+i, index);
                // maximum(0, lower_bound - value)
                l1_inf += fmax(0.0, tmp_bound-tmp);
                // printf("lower bounds: bound: %.4e, value: %.4e, result: %.4e\n", tmp_bound, tmp, fmax(0.0, tmp_bound-tmp));
            }
            // check upper bound
            mask_value = BLASFEO_DVECEL(qp_in->d_mask+i, nb[i] + ng[i] + index);
            if (mask_value == 1.0)
            {
                // upper bounds have the wrong sign! 
                // it is lower_bounds <= value <= -upper_bounds, therefore plus below
                tmp_bound = BLASFEO_DVECEL(qp_in->d+i, nb[i] + ng[i] + index);
                // printf("upper bounds: bound: %.4e, value: %.4e, result: %.4e\n", tmp_bound, tmp, fmax(0.0, tmp_bound-tmp));
                l1_inf += fmax(0.0, tmp_bound+tmp);
            }
        }
    }

    return l1_inf;
}

// /*
// This function calculates the l1 infeasibility by calculating the matrix vector product of the 
// constraints
// */
// static double manually_calculate_slacked_qp_l1_infeasibility(ocp_nlp_dims *dims, ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_sqp_wfqp_workspace *work, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
// {
//     int N = dims->N;
//     int *nx = dims->nx;
//     int *nu = dims->nu;
//     int *ns = dims->ns;
//     int *nns = mem->nns;
//     int *ni = dims->ni;
//     int *nh = dims->ni;

//     int *nb = qp_in->dim->nb;
//     int *ng = qp_in->dim->ng; //number of general two sided constraints

//     double l1_inf = 0.0;
//     int i, j;
//     double tmp, tmp_bound, mask_value;
//     int nominal_dims;

//     //TODO: Write this function more efficiently. Atm too many operations
//     for (i = 0; i <= N; i++)
//     {
//         // bounds on states and controls
//         for (j=0; j<nb[i]; ++j)
//         {
//             tmp = BLASFEO_DVECEL(qp_out->ux+i, qp_in->idxb[i][j]);
//             blasfeo_dvecin1(tmp, &work->nlp_work->tmp_ni+i, j);
//             blasfeo_dvecin1(tmp, &work->nlp_work->tmp_ni+i, nb[i]+ng[i]+j);
//         }
//         // general linear / linearized!
//         // tmp_ni = D * u + C * x
//         // lower bounds --> this seems to be correct and in accordance with slack variables
//         blasfeo_dgemv_t(nu[i]+nx[i], ng[i], 1.0, qp_in->DCt+i, 0, 0, qp_out->ux+i, 0,
//                         0.0, qp_in->d+i, nb[i], &work->nlp_work->tmp_ni+i, nb[i]);
//         blasfeo_dveccp(ng[i], &work->nlp_work->tmp_ni+i, nb[i], &work->nlp_work->tmp_ni+i, 2*nb[i]+ng[i]);
        
//         // upper bounds (seems to be correct but I do not understand why??)
//         // the sign of upper bound d is wrong!! We should use -d. Why is that?
//         // blasfeo_dgemv_t(nu[i]+nx[i], ng[i], 1.0, qp_in->DCt+i, 0, 0, qp_out->ux+i, 0,
//         //                 0.0, qp_in->d+i, 2*nb[i]+ng[i], &work->nlp_work->tmp_ni+i, 2*nb[i]+ng[i]);
//         for (j=0; j<2*nb[i]+2*ng[i]; ++j)
//         {
//             mask_value = BLASFEO_DVECEL(qp_in->d_mask+i, j);
//             if (mask_value == 1.0)
//             {
//                 tmp = BLASFEO_DVECEL(&work->nlp_work->tmp_ni+i, j);
//                 tmp_bound = BLASFEO_DVECEL(qp_in->d+i, j);
//                 if (j < nb[i] + ng[i])
//                 {
//                     // maximum(0, lower_bound - value)
//                     l1_inf += fmax(0.0, tmp_bound-tmp);
//                     // printf("lower bounds: bound: %.4e, value: %.4e, result: %.4e\n", tmp_bound, tmp, fmax(0.0, tmp_bound-tmp));
//                 }
//                 else
//                 {
//                     // upper bounds have the wrong sign! 
//                     // it is lower_bounds <= value <= -upper_bounds, therefore plus below
//                     // printf("upper bounds: value: %.4e, value: %.4e, result: %.4e\n", tmp_bound, tmp, fmax(0.0, tmp_bound+tmp));
//                     l1_inf += fmax(0.0, tmp_bound+tmp);
//                 }
//             }
//         }

//         // Printing
//         // printf("Manual multiplication at stage %d\n", i);
//         // printf("i=%d\n", i);
//         // printf("Print the lower bound vector\n");
//         // printf("product\n");
//         // blasfeo_print_dvec(ng[i], &work->nlp_work->tmp_ni+i, nb[i]);
//         // printf("bound\n");
//         // blasfeo_print_dvec(ng[i], qp_in->d+i, nb[i]);
//         // printf("d_mask\n");
//         // blasfeo_print_dvec(ng[i], qp_in->d_mask+i, nb[i]);
//         // printf("Print the upper bound vector\n");
//         // printf("product\n");
//         // blasfeo_print_dvec(ng[i], &work->nlp_work->tmp_ni+i, 2*nb[i]+ng[i]);
//         // printf("bound\n");
//         // blasfeo_print_dvec(ng[i], qp_in->d+i, 2*nb[i]+ng[i]);
//         // printf("d_mask\n");
//         // blasfeo_print_dvec(ng[i], qp_in->d_mask+i, 2*nb[i]+ng[i]);

//     }

//     return l1_inf;
// }


static void set_non_slacked_l1_penalties(ocp_nlp_config *config, ocp_nlp_dims *dims,
    ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_sqp_wfqp_memory *mem,
    ocp_nlp_workspace *work)
{
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int *nns = mem->nns;
    ocp_qp_in *qp_in = mem->nlp_mem->qp_in;

    // be aware of rqz_QP = [r, q, zl_NLP, zl_QP, zu_NLP, zu_QP]
    for (int stage = 0; stage <= dims->N; stage++)
    {
        // zl_QP
        blasfeo_dvecse(nns[stage], 1.0, qp_in->rqz+stage, nu[stage]+nx[stage]+ns[stage]);
        // zu_QP
        blasfeo_dvecse(nns[stage], 1.0, qp_in->rqz+stage, nu[stage]+nx[stage]+2*ns[stage]+nns[stage]);
        // printf("qp_in->rqz %d\n", stage);
        // blasfeo_print_exp_tran_dvec(nu[stage] +nx[stage] + 2*(nns[stage]+ns[stage]), qp_in->rqz+stage, 0);
    }
}


static void set_non_slacked_l2_penalties(ocp_nlp_config *config, ocp_nlp_dims *dims,
    ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_sqp_wfqp_memory *mem,
    ocp_nlp_workspace *work)
{
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

/*
Calculates the norm of the search direction, but the additional slack directions are removed.
If the problem is infeasibe, or the algorithm converges towards an infeasible point,
then the step size would not converge to 0 for our additional slack variables (since we do not do a delta 
update in the master problem).
*/
static double slacked_qp_out_compute_primal_nrm_inf(ocp_qp_out* qp_out, ocp_nlp_dims *dims, ocp_nlp_sqp_wfqp_memory* mem)
{
    double res = 0;
    double res_stage = 0;
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int *nns = mem->nns;

    for (int i = 0; i <= N; i++)
    {
        blasfeo_dvecnrm_inf(nx[i]+nu[i]+ns[i], qp_out->ux+i, 0, &res_stage);
        res += res_stage;
        blasfeo_dvecnrm_inf(ns[i], qp_out->ux+i, nx[i]+nu[i]+ns[i]+nns[i], &res_stage);
        res += res_stage;
    }
    return res;
}


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
    ocp_qp_out *qp_out = nlp_mem->qp_out;
    // solver_mem is not used in this function, but needed for DDP
    // the function is used in the config->globalization->step_update
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nz = dims->nz;
    int *ns = dims->ns;
    int *nns = mem->nns;
    int n_nominal_ineq_nlp, two_n_nominal_ineq_nlp;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // Assuming the variable order
        // [u x lb_slack ub_slack]
        // Assuming the variable order for QP direction
        // [u x lb_slack lb_relaxation_slack ub_slack ub_relaxation_slack]

        // step in primal variables
        blasfeo_daxpy(nx[i]+nu[i]+ns[i], alpha, qp_out->ux + i, 0, out_start->ux + i, 0, out_destination->ux + i, 0);
        blasfeo_daxpy(ns[i], alpha, qp_out->ux + i, nx[i]+nu[i]+ns[i]+nns[i],
                      out_start->ux + i, nx[i]+nu[i]+ns[i],
                      out_destination->ux + i, nx[i]+nu[i]+ns[i]);

        // update dual variables
        n_nominal_ineq_nlp = dims->ni[i] - ns[i];
        two_n_nominal_ineq_nlp = 2*n_nominal_ineq_nlp;
        // qp_out->lam = [lbu, ubu, lbx, ubx, lbg, ubg, lbh, ubh, lbs_NLP, lbs_QP, ubs_NLP, ubs_QP]
        // nlp_out->lam = [lbu, ubu, lbx, ubx, lbg, ubg, lbh, ubh, lbs_NLP, ----, ubs_NLP, ---]
        if (full_step_dual)
        {
            blasfeo_dveccp(two_n_nominal_ineq_nlp+ns[i], qp_out->lam+i, 0, out_destination->lam+i, 0);
            blasfeo_dveccp(ns[i], qp_out->lam+i, two_n_nominal_ineq_nlp+ns[i]+mem->nns[i],
                           out_destination->lam+i, two_n_nominal_ineq_nlp+ns[i]);
            if (i < N)
            {
                blasfeo_dveccp(nx[i+1], qp_out->pi+i, 0, out_destination->pi+i, 0);
            }
        }
        else
        {
            // update duals with alpha step
            blasfeo_daxpby(two_n_nominal_ineq_nlp+ns[i], 1.0-alpha, out_start->lam+i, 0, alpha, qp_out->lam+i, 0, out_destination->lam+i, 0);
            blasfeo_daxpby(ns[i], 1.0-alpha, out_start->lam+i, two_n_nominal_ineq_nlp+ns[i]+mem->nns[i], alpha,
                           qp_out->lam+i, two_n_nominal_ineq_nlp+ns[i]+mem->nns[i], out_destination->lam+i, two_n_nominal_ineq_nlp+ns[i]);

            if (i < N)
            {
                blasfeo_daxpby(nx[i+1], 1.0-alpha, out_start->pi+i, 0, alpha, qp_out->pi+i, 0, out_destination->pi+i, 0);
            }
        }

        // linear update of algebraic variables using state and input sensitivity
        if (i < N)
        {
            // out->z = nlp_mem->z_alg + alpha * dzdux * qp_out->ux
            blasfeo_dgemv_t(nu[i]+nx[i], nz[i], alpha, nlp_mem->dzduxt+i, 0, 0,
                    qp_out->ux+i, 0, 1.0, nlp_mem->z_alg+i, 0, out_destination->z+i, 0);
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

/************************************************
 * functions for QP preparation
 ************************************************/

static void ocp_nlp_sqp_wfqp_prepare_hessian_evaluation(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
    ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_workspace *work)
{
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    int N = dims->N;

    // int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;

    if (!nlp_mem->compute_hess)
    {
        printf("ocp_nlp_sqp_wfqp_prepare_hessian_evaluation: constant hessian not supported!\n\n");
        exit(1);
    }

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // TODO: first compute cost hessian (without adding) and avoid setting everything to zero?
        // init Hessians to 0

        // TODO: avoid setting qp_in->RSQ to zero in ocp_nlp_approximate_qp_matrices?
        blasfeo_dgese(nu[i] + nx[i], nu[i] + nx[i], 0.0, mem->RSQ_constr+i, 0, 0);
        blasfeo_dgese(nu[i] + nx[i], nu[i] + nx[i], 0.0, mem->RSQ_cost+i, 0, 0);
    }

    for (int i = 0; i < N; i++)
    {
        config->dynamics[i]->memory_set_RSQrq_ptr(mem->RSQ_constr+i, nlp_mem->dynamics[i]);
    }
    for (int i = 0; i <= N; i++)
    {
        config->cost[i]->memory_set_RSQrq_ptr(mem->RSQ_cost+i, nlp_mem->cost[i]);
        config->constraints[i]->memory_set_RSQrq_ptr(mem->RSQ_constr+i, nlp_mem->constraints[i]);
    }
    return;
}

// update QP rhs for SQP (step prim var, abs dual var)
// - use cost gradient and dynamics residual from memory
// - evaluate constraints wrt bounds -> allows to update all bounds between preparation and feedback phase.
void ocp_nlp_sqp_wfqp_approximate_qp_constraint_vectors(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
    ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_workspace *work)
{
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    int N = dims->N;

    // int *nv = dims->nv;
    int *nx = dims->nx;
    int *ns = dims->ns;
    int *ni = dims->ni;
    int *nns = mem->nns;


#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // b
        if (i < N)
            blasfeo_dveccp(nx[i + 1], nlp_mem->dyn_fun + i, 0, nlp_mem->qp_in->b + i, 0);

        // evaluate constraint residuals
        config->constraints[i]->update_qp_vectors(config->constraints[i], dims->constraints[i],
            in->constraints[i], opts->constraints[i], nlp_mem->constraints[i], work->constraints[i]);

        // copy ineq function value into nlp mem, then into QP
        struct blasfeo_dvec *ineq_fun = config->constraints[i]->memory_get_fun_ptr(nlp_mem->constraints[i]);
        blasfeo_dveccp(2 * ni[i], ineq_fun, 0, nlp_mem->ineq_fun + i, 0);

        // d
        int n_nominal_ineq_nlp = ni[i] - ns[i];

        // blasfeo_dveccp(2 * ni[i], nlp_mem->ineq_fun + i, 0, nlp_mem->qp_in->d + i, 0);
        blasfeo_dveccp(2*n_nominal_ineq_nlp+ns[i], nlp_mem->ineq_fun + i, 0, nlp_mem->qp_in->d + i, 0);
        blasfeo_dveccp(ns[i], nlp_mem->ineq_fun + i, 2*n_nominal_ineq_nlp+ns[i], nlp_mem->qp_in->d + i, 2*n_nominal_ineq_nlp+ns[i]+nns[i]);
    }
}


static void ocp_nlp_sqp_wfqp_setup_qp_objective(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
    ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_workspace *work, double objective_multiplier)
{
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    int N = dims->N;

    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int *nns = mem->nns;

    int nxu;
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        /* Hessian matrices */
        // hess_QP = objective_multiplier * hess_cost + hess_constraints
        nxu = nx[i]+nu[i];
        // TODO: axpby for matrices? Do we need to take slacks here into account as well? I.e., scale slack Hessian?
        if (objective_multiplier == 0.0)
        {
            // Either we use the exact objective Hessian
            // blasfeo_dgecp(nxu, nxu, mem->RSQ_constr+i, 0, 0, nlp_mem->qp_in->RSQrq+i, 0, 0);
            // blasfeo_dgead(nxu, nxu, objective_multiplier, mem->RSQ_cost+i, 0, 0, nlp_mem->qp_in->RSQrq+i, 0, 0); I think we do not need this here
            
            // We use the identity matrix Hessian
            blasfeo_dgese(nxu, nxu, 0.0, nlp_mem->qp_in->RSQrq+i, 0, 0);
            blasfeo_ddiare(nxu, 1e-4, nlp_mem->qp_in->RSQrq+i, 0, 0);  // dPsi_dx is unit now
        }
        else
        {
            blasfeo_dgecp(nxu, nxu, mem->RSQ_constr+i, 0, 0, nlp_mem->qp_in->RSQrq+i, 0, 0);
            //
            blasfeo_dgead(nxu, nxu, objective_multiplier, mem->RSQ_cost+i, 0, 0, nlp_mem->qp_in->RSQrq+i, 0, 0);
        }
        // Z -- slack matrix
        blasfeo_dveccpsc(ns[i], objective_multiplier, mem->Z_cost_module+i, 0, nlp_mem->qp_in->Z+i, 0);
        blasfeo_dveccpsc(ns[i], objective_multiplier, mem->Z_cost_module+i, 0, nlp_mem->qp_in->Z+i, ns[i]+mem->nns[i]);

        /* vectors */
        // g
        blasfeo_dveccpsc(nx[i]+nu[i]+ns[i], objective_multiplier, nlp_mem->cost_grad + i, 0, nlp_mem->qp_in->rqz + i, 0);
        blasfeo_dveccpsc(ns[i], objective_multiplier, nlp_mem->cost_grad + i, nx[i]+nu[i]+ns[i], nlp_mem->qp_in->rqz + i, nx[i]+nu[i]+ns[i]+nns[i]);
    }
}

/*
Solves the QP. We either solve the feasibility QP or the standard l1-relaxed QP
*/
static int prepare_and_solve_QP(ocp_nlp_config* config, ocp_nlp_sqp_wfqp_opts* opts, ocp_qp_in* qp_in, ocp_qp_out* qp_out,
                    ocp_nlp_dims *dims, ocp_nlp_sqp_wfqp_memory* mem, ocp_nlp_in* nlp_in, ocp_nlp_out* nlp_out,
                    ocp_nlp_memory* nlp_mem, ocp_nlp_workspace* nlp_work, int sqp_iter, bool solve_feasibility_qp,
                    acados_timer timer0, acados_timer timer1)
{
    ocp_nlp_opts* nlp_opts = opts->nlp_opts;
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    // ocp_nlp_res *nlp_res = nlp_mem->nlp_res;
    ocp_nlp_timings *nlp_timings = nlp_mem->nlp_timings;

    // (typically) no warm start at first iteration
    if (sqp_iter == 0 && !opts->warm_start_first_qp)
    {
        int tmp_int = 0;
        qp_solver->opts_set(qp_solver, nlp_opts->qp_solver_opts, "warm_start", &tmp_int);
    }
    // Load input to QP and regularize Hessian
    if (solve_feasibility_qp)
    {
        ocp_nlp_sqp_wfqp_setup_qp_objective(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work, 0.0);
    }
    else
    {
        ocp_nlp_sqp_wfqp_setup_qp_objective(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work, nlp_mem->objective_multiplier);
    }
    // TODO: if we solve the feasibility QP, we probably do not need or want the LM term?
    ocp_nlp_add_levenberg_marquardt_term(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, mem->alpha, sqp_iter);
    // regularize Hessian
    acados_tic(&timer1);
    config->regularize->regularize(config->regularize, dims->regularize,
                                            nlp_opts->regularize, nlp_mem->regularize_mem);
    nlp_timings->time_reg += acados_toc(&timer1);
    // Show input to QP
    if (nlp_opts->print_level > 3)
    {
        printf("\n\nSQP: ocp_qp_in at iteration %d\n", sqp_iter);
        print_ocp_qp_in(qp_in);
    }

#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
    ocp_nlp_dump_qp_in_to_file(qp_in, sqp_iter, 0);
#endif

    int qp_status = ocp_nlp_solve_qp_and_correct_dual(config, dims, nlp_opts, nlp_mem, nlp_work, false, NULL, qp_out);

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
    int qp_iter = qp_info_->num_iter;

    // save statistics of last qp solver call
    if (sqp_iter+1 < mem->stat_m && !solve_feasibility_qp)
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
    // TODO: @david: How to log statistics?
    // We solve 2 QPs -> need 2 times qp_iter and qp_status?
    // Or do we make one row for one QP solve?
    return qp_status;
}

/************************************************
 * residual functions
 ************************************************/

void ocp_nlp_sqp_wfqp_res_compute(ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_res *res,
                         ocp_nlp_sqp_wfqp_memory *mem)
{
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    // extract dims
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;

    double tmp_res;
    double tmp;

    // res_stat
    for (int i = 0; i <= N; i++)
    {
        // blasfeo_daxpy(nv[i], -nlp_mem->objective_multiplier, nlp_mem->cost_grad + i, 0, nlp_mem->ineq_adj + i, 0,
        //               res->res_stat + i, 0);
        blasfeo_daxpy(nv[i], -1.0, nlp_mem->cost_grad + i, 0, nlp_mem->ineq_adj + i, 0,
                      res->res_stat + i, 0);
        blasfeo_daxpy(nu[i] + nx[i], 1.0, nlp_mem->dyn_adj + i, 0, res->res_stat + i, 0,
                      res->res_stat + i, 0);
        blasfeo_dvecnrm_inf(nv[i], res->res_stat + i, 0, &tmp_res);
        blasfeo_dvecse(1, tmp_res, &res->tmp, i);
    }
    blasfeo_dvecnrm_inf(N+1, &res->tmp, 0, &res->inf_norm_res_stat);

    // res_eq
    for (int i = 0; i < N; i++)
    {
        blasfeo_dveccp(nx[i + 1], nlp_mem->dyn_fun + i, 0, res->res_eq + i, 0);
        blasfeo_dvecnrm_inf(nx[i + 1], res->res_eq + i, 0, &tmp_res);
        blasfeo_dvecse(1, tmp_res, &res->tmp, i);
    }
    blasfeo_dvecnrm_inf(N, &res->tmp, 0, &res->inf_norm_res_eq);

    // res_ineq
    res->inf_norm_res_ineq = 0.0;
    for (int i = 0; i <= N; i++)
    {
        for (int j=0; j<2*ni[i]; j++)
        {
            tmp = BLASFEO_DVECEL(nlp_mem->ineq_fun+i, j);
            if (tmp > res->inf_norm_res_ineq)
            {
                res->inf_norm_res_ineq = tmp;
            }
        }
    }

    // res_comp
    res->inf_norm_res_comp = 0.0;
    for (int i = 0; i <= N; i++)
    {
        blasfeo_dvecmul(2 * ni[i], out->lam + i, 0, nlp_mem->ineq_fun+i, 0, res->res_comp + i, 0);
        blasfeo_dvecnrm_inf(2 * ni[i], res->res_comp + i, 0, &tmp_res);
        blasfeo_dvecse(1, tmp_res, &res->tmp, i);
    }
    blasfeo_dvecnrm_inf(N+1, &res->tmp, 0, &res->inf_norm_res_comp);
}

/************************************************
* MAIN OPTIMIZATION ROUTINE
************************************************/
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
    // ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
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
    mem->nlp_mem->status = ACADOS_READY;
    nlp_mem->objective_multiplier = opts->initial_objective_multiplier;


#if defined(ACADOS_WITH_OPENMP)
    // backup number of threads
    int num_threads_bkp = omp_get_num_threads();
    // set number of threads
    omp_set_num_threads(opts->nlp_opts->num_threads);
#endif

    ocp_nlp_initialize_submodules(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
    set_non_slacked_l2_penalties(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work);
    set_non_slacked_l1_penalties(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work);

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
            ocp_nlp_sqp_wfqp_prepare_hessian_evaluation(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work);
            acados_tic(&timer1);
            ocp_nlp_approximate_qp_matrices(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
            ocp_nlp_sqp_wfqp_approximate_qp_constraint_vectors(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work);

            if (nlp_opts->with_adaptive_levenberg_marquardt || config->globalization->needs_objective_value() == 1)
            {
                ocp_nlp_get_cost_value_from_submodules(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
            }
            //

            nlp_timings->time_lin += acados_toc(&timer1);

            // compute nlp residuals; TODO: double check
            // ocp_nlp_sqp_wfqp_res_compute(dims, nlp_in, nlp_out, nlp_res, mem);
            ocp_nlp_res_compute(dims, nlp_in, nlp_out, nlp_res, nlp_mem);
            ocp_nlp_res_get_inf_norm(nlp_res, &nlp_out->inf_norm_res);
        }

        // TODO: @david: where do we need this? move somewhere else?
        set_slack_variable_values(config, dims, mem);
        // printf("Current slack values of non-original NLP slacks:\n");
        // for (int jj=0; jj<=dims->N;++jj)
        // {
        //     blasfeo_print_dvec(mem->nns[jj], mem->slacks_not_in_original_nlp + jj, 0);
        // }

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

        // printf("Current iterate: \n");
    	// for (int jj=0; jj<=dims->N;++jj)
        // {
        //     blasfeo_print_dvec(dims->nx[jj]+dims->nu[jj], nlp_out->ux + jj, 0);
        // }

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
        
        double multiplier_norm_inf = get_multiplier_norm_inf(dims, qp_out);
        printf("Multiplier norm inf is: %.4e\n", multiplier_norm_inf);

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

        /* solve 1. QP: We solve without gradient and only with constraint Hessian */
        qp_status = prepare_and_solve_QP(config, opts, qp_in, nlp_work->tmp_qp_out, dims, mem, nlp_in, nlp_out,
                    nlp_mem, nlp_work, sqp_iter, true, timer0, timer1);
        if (qp_status != ACADOS_SUCCESS)
        {
            if (nlp_opts->print_level > 1)
            {
                printf("\nFailure in QP 1, got status %d!\n", qp_status);
            }
            nlp_mem->status = qp_status;
            nlp_mem->iter = sqp_iter;
            nlp_timings->time_tot = acados_toc(&timer0);
#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif
            return nlp_mem->status;
        }

        /* solve 2. QP: We solve the standard l1-relaxed QP with gradient */
        qp_status = prepare_and_solve_QP(config, opts, qp_in, qp_out, dims, mem, nlp_in, nlp_out,
                    nlp_mem, nlp_work, sqp_iter, false, timer0, timer1);
        if (qp_status != ACADOS_SUCCESS)
        {
            if (nlp_opts->print_level > 1)
            {
                printf("\nFailure in QP 2, got status %d!\n", qp_status);
            }
            nlp_mem->status = qp_status;
            nlp_mem->iter = sqp_iter;
            nlp_timings->time_tot = acados_toc(&timer0);
#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif
            return nlp_mem->status;
        }
        //TODO: We should some day also account for, if the QP was not solved, e.g., max iter was reached

        // @david: now you have two directions:
        // nlp_work->tmp_qp_out = d without cost;
        // nlp_mem->qp_out = d with objective_multiplier

        multiplier_norm_inf = get_multiplier_norm_inf(dims, nlp_work->tmp_qp_out);
        printf("Feasibility Multiplier norm inf is: %.4e\n", multiplier_norm_inf);
        multiplier_norm_inf = get_multiplier_norm_inf(dims, qp_out);
        printf("Optimality Multiplier norm inf is: %.4e\n", multiplier_norm_inf);

        double current_l1_infeasibility = ocp_nlp_get_l1_infeasibility(config, dims, nlp_mem);
        printf("Current l1 infeasibility: %.4e\n", current_l1_infeasibility);

        // Calculate linearized l1-infeasibility for d_steering
        // double l1_inf_QP_feasibility = get_slacked_qp_l1_infeasibility(dims, mem, nlp_work->tmp_qp_out);
        // printf("linearized l1_inf_feas: %.4e\n", l1_inf_QP_feasibility);
        double manual_l1_inf_QP_feasibility = manually_calculate_slacked_qp_l1_infeasibility(dims, mem, work, nlp_mem->qp_in, nlp_work->tmp_qp_out);
        printf("manual l1_inf_feas: %.4e\n", manual_l1_inf_QP_feasibility);

        // Calculate linearized l1-infeasibility for d_predictor
        // double l1_inf_QP_optimality = get_slacked_qp_l1_infeasibility(dims, mem, nlp_mem->qp_out);
        // printf("linearized l1_inf_opt: %.4e\n", l1_inf_QP_optimality);
        double manual_l1_inf_QP_optimality = manually_calculate_slacked_qp_l1_infeasibility(dims, mem, work, nlp_mem->qp_in, nlp_mem->qp_out);
        printf("manual l1_inf_opt: %.4e\n", manual_l1_inf_QP_optimality);

        // predicted infeasibility reduction of feasibility QP should always be non-negative
        double pred_l1_inf_QP_feasibility, pred_l1_inf_QP_optimality;
        if (current_l1_infeasibility < fmin(opts->tol_ineq, opts->tol_eq)) // do this to avoid some weird negative values in pred
        {
            pred_l1_inf_QP_feasibility = 0.0;
            pred_l1_inf_QP_optimality = 0.0;
        }
        else
        {
            pred_l1_inf_QP_feasibility = current_l1_infeasibility - manual_l1_inf_QP_feasibility;
            pred_l1_inf_QP_optimality = current_l1_infeasibility - manual_l1_inf_QP_optimality;
        }
        pred_l1_inf_QP_feasibility = calculate_predicted_l1_inf_reduction(opts, current_l1_infeasibility, manual_l1_inf_QP_feasibility);
        pred_l1_inf_QP_optimality = calculate_predicted_l1_inf_reduction(opts, current_l1_infeasibility, manual_l1_inf_QP_optimality);
        printf("pred_l1_inf_QP_feasibility: %.4e\n", pred_l1_inf_QP_feasibility);
        printf("pred_l1_inf_QP_optimality: %.4e\n", pred_l1_inf_QP_optimality);

        // It seems appropriate that the fraction for sufficient improvement
        // in infeasibility is adaptive. So, if inf is large, and the improvement is small
        // relative to infeasibility we should have a direction that is closer to 
        // the feasibility direction??
        double kappa = calculate_search_direction_interpolation_factor(opts,
                                                                pred_l1_inf_QP_feasibility,
                                                                pred_l1_inf_QP_optimality,
                                                                manual_l1_inf_QP_optimality,
                                                                manual_l1_inf_QP_feasibility);

        // Calculate search direction
        setup_search_direction(mem, dims, qp_out, nlp_work->tmp_qp_out, qp_out, kappa);
        
        double pred_l1_inf_search_direction = calculate_predicted_l1_inf_reduction(opts, current_l1_infeasibility, manually_calculate_slacked_qp_l1_infeasibility(dims, mem, work, nlp_mem->qp_in, nlp_mem->qp_out));
        printf("pred_l1_inf_search_direction: %.4e\n", pred_l1_inf_search_direction);
        //---------------------------------------------------------------------
        // scale_multiplier(dims, nlp_mem, qp_out);

        // Calculate optimal QP objective (needed for globalization)
        if (config->globalization->needs_qp_objective_value() == 1)
        {
            nlp_mem->qp_cost_value = ocp_nlp_compute_qp_objective_value(dims, qp_in, qp_out, nlp_work);
            // is this correct?
            if (kappa == 1.0)
            {
                nlp_mem->predicted_infeasibility_reduction = pred_l1_inf_QP_optimality;
            }
            else
            {
                nlp_mem->predicted_infeasibility_reduction = pred_l1_inf_search_direction;
            }
        }

        // Compute the step norm
        if (opts->tol_min_step_norm > 0.0 || nlp_opts->log_primal_step_norm)
        {
            // For the moment we do not care about the artificial slack variables to keep problem
            // feasible
            mem->step_norm = slacked_qp_out_compute_primal_nrm_inf(qp_out, dims, mem);
            // printf("Step norm: %.4e\n", mem->step_norm);
            if (nlp_opts->log_primal_step_norm)
                mem->primal_step_norm[sqp_iter] = mem->step_norm;
        }
        /* end solve QP */
        if (current_l1_infeasibility > fmin(opts->tol_eq, opts->tol_ineq) && mem->step_norm < opts->tol_min_step_norm)
        {
            printf("Problems seems to be converged to an infeasible stationary point!\n");
            nlp_mem->status = ACADOS_INFEASIBLE;
            return nlp_mem->status;
        }

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
            nlp_mem->status = globalization_status;
            nlp_mem->iter = sqp_iter;
            nlp_timings->time_tot = acados_toc(&timer0);
#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif
            return nlp_mem->status;
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
    return nlp_mem->status;
}


// TODO(david): use common?
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
        nns = mem->nns[stage];
        for (int i=0; i<nns; i++)
        {
            idxs_rev[idxns[i]] = i+dims->ns[stage];
        }
        // printf("set qp_in->idxs_rev to\n");
        // for (int i=0; i<dims->ni[stage]-dims->ns[stage]; i++)
        //     printf("%d\t", idxs_rev[i]);
        // printf("\n");
    }

    ocp_nlp_precompute_common(config, dims, nlp_in, nlp_out, opts->nlp_opts, nlp_mem, nlp_work);

    // overwrite output pointers normally set in ocp_nlp_alias_memory_to_submodules
    for (int stage = 0; stage <= dims->N; stage++)
    {
        config->cost[stage]->memory_set_Z_ptr(mem->Z_cost_module+stage, nlp_mem->cost[stage]);
    }

    return ACADOS_SUCCESS;
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
