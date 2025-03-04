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
    opts->eval_residual_at_max_iter = true; // we want to know in last iteration if we converged
    opts->use_QP_l1_inf_from_slacks = false; // if manual calculation used, results seem more accurate and solver performs better!

    opts->use_exact_hessian_in_feas_qp = false;
    opts->search_direction_mode = NOMINAL_QP;
    opts->watchdog_zero_slacks_max = 2;
    opts->allow_direction_mode_switch = true;

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
        else if (!strcmp(field, "use_exact_hessian_in_feas_qp"))
        {
            bool* use_exact_hessian_in_feas_qp = (bool *) value;
            opts->use_exact_hessian_in_feas_qp = *use_exact_hessian_in_feas_qp;
        }
        else if (!strcmp(field, "search_direction_mode"))
        {
            bool* search_direction_mode = (bool *) value;
            opts->search_direction_mode = *search_direction_mode;
        }
        else if (!strcmp(field, "allow_direction_mode_switch"))
        {
            bool* allow_direction_mode_switch = (bool *) value;
            opts->allow_direction_mode_switch = *allow_direction_mode_switch;
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

acados_size_t ocp_nlp_sqp_wfqp_memory_calculate_size(void *config_, void *dims_, void *opts_, void *in_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_in *in = in_;
    ocp_nlp_sqp_wfqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    acados_size_t size = 0;
    int N = dims->N;

    size += sizeof(ocp_nlp_sqp_wfqp_memory);

    // nlp mem
    size += ocp_nlp_memory_calculate_size(config, dims, nlp_opts, in);

    // relaxed QP solver
    size += ocp_qp_xcond_solver_memory_calculate_size(config->relaxed_qp_solver, dims->relaxed_qp_solver, opts->nlp_opts->qp_solver_opts);
    size += ocp_qp_xcond_solver_workspace_calculate_size(config->relaxed_qp_solver, dims->relaxed_qp_solver, opts->nlp_opts->qp_solver_opts);
    size += ocp_qp_in_calculate_size(dims->relaxed_qp_solver->orig_dims);
    size += ocp_qp_out_calculate_size(dims->relaxed_qp_solver->orig_dims);

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
    // nlp_idxs_rev
    size += (N + 1) * sizeof(int *);

    int nns, nsbu, nbu, nsbx, nbx, n_nominal_ineq_nlp, stage;
    for (stage = 0; stage <= N; stage++)
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
        // nlp_idxs_rev
        size += (dims->nb[stage] + dims->ng[stage] + dims->ni_nl[stage]) * sizeof(int);

        // multipliers for the feasibility QP
        size += 1 * blasfeo_memsize_dvec(2 * dims->ni[stage]);  // lam_feasibility
        if (stage < N)
        {
            size += 1 * blasfeo_memsize_dvec(dims->nx[stage + 1]);  // pi_feasibility
        }
        size += 1 * blasfeo_memsize_dvec(dims->nv[stage]);      // res_stat_feasibility

        // Z_cost_module
        size += blasfeo_memsize_dvec(2*dims->ns[stage]);

        // RSQ_cost, RSQ_constr
        size += 2*blasfeo_memsize_dmat(dims->nx[stage]+dims->nu[stage], dims->nx[stage]+dims->nu[stage]);
    }
    // nns
    size += (N+1) * sizeof(int);
    // multipliers for the feasibility QP
    size += 1 * (N + 1) * sizeof(struct blasfeo_dvec);  // lam_feasibility
    size += 1 * N * sizeof(struct blasfeo_dvec);  // pi_feasibility
    size += 1 * (N + 1) * sizeof(struct blasfeo_dvec);  // res_stat_feasibility
    // Z_cost_module
    size += (N + 1) * sizeof(struct blasfeo_dvec);
    // RSQ_cost, RSQ_constr
    size += 2*(N + 1) * sizeof(struct blasfeo_dmat);

    size += 3*8;  // align
    size += 64;  // blasfeo_mem align

    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_nlp_sqp_wfqp_memory_assign(void *config_, void *dims_, void *opts_, void *in_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_in *in = in_;
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
    mem->nlp_mem = ocp_nlp_memory_assign(config, dims, nlp_opts, in, c_ptr);
    c_ptr += ocp_nlp_memory_calculate_size(config, dims, nlp_opts, in);

    /* relaxed QP solver */
    // mem
    mem->relaxed_qp_solver_mem = ocp_qp_xcond_solver_memory_assign(config->relaxed_qp_solver, dims->relaxed_qp_solver, opts->nlp_opts->qp_solver_opts, c_ptr);
    c_ptr += ocp_qp_xcond_solver_memory_calculate_size(config->relaxed_qp_solver, dims->relaxed_qp_solver, opts->nlp_opts->qp_solver_opts);
    // work
    mem->relaxed_qp_solver_work = (ocp_qp_xcond_solver_workspace*) c_ptr;
    c_ptr += ocp_qp_xcond_solver_workspace_calculate_size(config->relaxed_qp_solver, dims->relaxed_qp_solver, opts->nlp_opts->qp_solver_opts);
    // in & out
    mem->relaxed_qp_in = ocp_qp_in_assign(dims->relaxed_qp_solver->orig_dims, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(dims->relaxed_qp_solver->orig_dims);
    mem->relaxed_qp_out = ocp_qp_out_assign(dims->relaxed_qp_solver->orig_dims, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(dims->relaxed_qp_solver->orig_dims);
    // nominal
    mem->relaxed_qp_solver.config = config->relaxed_qp_solver;
    mem->relaxed_qp_solver.dims = dims->relaxed_qp_solver;
    mem->relaxed_qp_solver.opts = opts->nlp_opts->qp_solver_opts;
    mem->relaxed_qp_solver.mem = mem->relaxed_qp_solver_mem;
    mem->relaxed_qp_solver.work = mem->relaxed_qp_solver_work;

    // pi_feasibility
    assign_and_advance_blasfeo_dvec_structs(N, &mem->pi_feasibility, &c_ptr);
    // lam_feasibility
    assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->lam_feasibility, &c_ptr);
    // res_stat_feasibility
    assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->res_stat_feasibility, &c_ptr);

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

    // assign_and_advance_int_ptrs(N+1, &mem->nlp_idxs_rev, &c_ptr);

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

        // assign_and_advance_int(dims->nb[stage]+dims->ng[stage]+dims->ni_nl[stage], &mem->nlp_idxs_rev[stage], &c_ptr);
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
    // pi_feasibility
    for (int i = 0; i < N; ++i)
    {
        assign_and_advance_blasfeo_dvec_mem(dims->nx[i + 1], mem->pi_feasibility + i, &c_ptr);
    }
    // lam_feasibility
    for (int i = 0; i <= N; ++i)
    {
        assign_and_advance_blasfeo_dvec_mem(2 * dims->ni[i], mem->lam_feasibility + i, &c_ptr);
    }
    // res_stat_feasibility
    for (int i = 0; i <= N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(dims->nv[i], mem->res_stat_feasibility + i, &c_ptr);
    }
    // Z_cost_module
    for (int i = 0; i <= N; ++i)
    {
        assign_and_advance_blasfeo_dvec_mem(2*dims->ns[i], mem->Z_cost_module + i, &c_ptr);
    }
    assert((char *) raw_memory + ocp_nlp_sqp_wfqp_memory_calculate_size(config, dims, opts, in) >= c_ptr);

    // initialize with zero
    for(int i=0; i<N; i++)
    {
        blasfeo_dvecse(dims->nx[i+1], 0.0, mem->pi_feasibility+i, 0);
        blasfeo_dvecse(2*dims->ni[i], 0.0, mem->lam_feasibility+i, 0);
    }
    blasfeo_dvecse(2*dims->ni[N], 0.0, mem->lam_feasibility+N, 0);

    return mem;
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_nlp_sqp_wfqp_workspace_calculate_size(void *config_, void *dims_, void *opts_, void *in_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_in *in = in_;
    ocp_nlp_sqp_wfqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    acados_size_t size = 0;

    // sqp
    size += sizeof(ocp_nlp_sqp_wfqp_workspace);

    // nlp
    size += ocp_nlp_workspace_calculate_size(config, dims, nlp_opts, in);

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
         ocp_nlp_sqp_wfqp_opts *opts, ocp_nlp_in *in, ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_sqp_wfqp_workspace *work)
{
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    // sqp
    char *c_ptr = (char *) work;
    c_ptr += sizeof(ocp_nlp_sqp_wfqp_workspace);

    // nlp
    work->nlp_work = ocp_nlp_workspace_assign(config, dims, nlp_opts, in, nlp_mem, c_ptr);
    c_ptr += ocp_nlp_workspace_calculate_size(config, dims, nlp_opts, in);

    if (nlp_opts->ext_qp_res)
    {
        // qp res
        work->nlp_work->qp_res = ocp_qp_res_assign(dims->qp_solver->orig_dims, c_ptr);
        c_ptr += ocp_qp_res_calculate_size(dims->qp_solver->orig_dims);

        // qp res ws
        work->nlp_work->qp_res_ws = ocp_qp_res_workspace_assign(dims->qp_solver->orig_dims, c_ptr);
        c_ptr += ocp_qp_res_workspace_calculate_size(dims->qp_solver->orig_dims);
    }

    assert((char *) work + ocp_nlp_sqp_wfqp_workspace_calculate_size(config, dims, opts, in) >= c_ptr);

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
            printf("Stopped: Maximum iterations Reached.\n");
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
            printf("Optimal solution found! Converged to KKT point.\n");
        }
        return true;
    }

    // check for infeasible problem
    // if (nlp_res->inf_norm_res_eq < opts->tol_eq
    //     && nlp_res->inf_norm_res_ineq > opts->tol_ineq
    //     && mem->inf_norm_res_stat_feasibility < opts->tol_stat
    //     && mem->inf_norm_res_comp_feasibility < opts->tol_comp)
    // {
    //     mem->nlp_mem->status = ACADOS_INFEASIBLE;
    //     if (opts->nlp_opts->print_level > 0)
    //     {
    //         printf("Converged to infeasible stationary point! Problem might be locally infeasible!\n");
    //     }
    //     return true;
    // }

    // check for small step
    if (opts->tol_min_step_norm > 0.0 && (n_iter > 0) && (mem->step_norm < opts->tol_min_step_norm))
    {
        if (opts->nlp_opts->print_level > 0)
        {
            if (nlp_res->inf_norm_res_eq < opts->tol_eq && nlp_res->inf_norm_res_ineq < opts->tol_ineq)
            {
                printf("Stopped: Converged to feasible Point. Step size is < tol_eq.\n");
            }
            else
            {
                printf("Stopped: Converged to infeasible Point. Step size is < tol_eq.\n");
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
            printf("Stopped: Maximum iterations reached.\n");
        }
        return true;
    }

    // convergence to FJ point?

    return false;
}

/************************************************
 * output
 ************************************************/
static void print_iteration(int iter, ocp_nlp_config *config, ocp_nlp_res *nlp_res, ocp_nlp_sqp_wfqp_memory *mem,
    ocp_nlp_opts *nlp_opts, double prev_levenberg_marquardt, int qp_status, int qp_iter)
{
ocp_nlp_memory *nlp_mem = mem->nlp_mem;
// print iteration header
if (iter % 10 == 0)
{
ocp_nlp_common_print_iteration_header();
printf("%9s   %9s   %8s   ", "step_norm", "step_type", "lm_reg.");
config->globalization->print_iteration_header();
printf("\n");
}
// print iteration
ocp_nlp_common_print_iteration(iter, nlp_res);
printf("%9.2e   %9s   %8.2e   ", mem->step_norm, mem->search_direction_type, prev_levenberg_marquardt);
config->globalization->print_iteration(nlp_mem->cost_value, nlp_opts->globalization, nlp_mem->globalization);
printf("\n");
}

/************************************************
 * functions
 ************************************************/
/*
Compute infinity norm of QP multipliers in NLP space.
Assumes that the masked multipliers are always zero.
*/
static void compute_qp_multiplier_norm_inf(ocp_nlp_sqp_wfqp_memory* mem, ocp_nlp_dims *dims,
                                           ocp_qp_out *qp_out, ocp_qp_in *qp_in,
                                           bool get_feasibility_multipliers)
{
    int i,j;
    int N = dims->N;
    int *nx = dims->nx;
    double norm_pi = 0.0;
    double norm_lam_slacked_constraints = 0.0;
    double norm_lam_hard_constr = 0.0;

    // compute inf norm of pi
    for (i = 0; i < N; i++)
    {
        for (j=0; j<nx[i+1]; j++)
        {
            norm_pi = fmax(norm_pi, fabs(BLASFEO_DVECEL(qp_out->pi+i, j)));
        }
    }
    if (get_feasibility_multipliers)
    {
        mem->norm_feas_qp_pi = norm_pi;
    }
    else
    {
        mem->norm_opt_qp_pi = norm_pi;
    }

    /* structure of QP and NLP iterates: */
    // qp_out->lam = [lbu, lbx, lg, l_nl, ubu, ubx, ug, u_nl, lbs_NLP, lbs_QP, ubs_NLP, ubs_QP]
    // nlp_out->lam = [lbu, lbx, lg, l_nl, ubu, ubx, ug, u_nl, lbs_NLP, ----, ubs_NLP, ---]
    int n_unslacked_bounds, n_nominal_ineq_nlp;
    for (i = 0; i <= N; i++)
    {
        n_nominal_ineq_nlp = dims->nb[i]+dims->ng[i]+dims->ni_nl[i];
        if (i == 0)
        {
            // we do not slack the initial state conditions!
            n_unslacked_bounds = qp_in->dim->nbu[i] + qp_in->dim->nbx[i];
        }
        else
        {
            n_unslacked_bounds = qp_in->dim->nbu[i];
        }

        // Unslacked bounds multipliers
        for (j=0; j<n_unslacked_bounds; j++)
        {
            norm_lam_hard_constr = fmax(norm_lam_hard_constr, BLASFEO_DVECEL(qp_out->lam+i, j));
            norm_lam_hard_constr = fmax(norm_lam_hard_constr, BLASFEO_DVECEL(qp_out->lam+i, n_nominal_ineq_nlp+j));
        }

        for (j=n_unslacked_bounds; j<n_nominal_ineq_nlp; j++)
        {
            // lb_not_l1_slacked, lg, l_nl
            norm_lam_slacked_constraints = fmax(norm_lam_slacked_constraints, BLASFEO_DVECEL(qp_out->lam+i, j));
            // ub_not_l1_slacked, ug, u_nl
            norm_lam_slacked_constraints = fmax(norm_lam_slacked_constraints, BLASFEO_DVECEL(qp_out->lam+i, j+n_nominal_ineq_nlp));
        }
        // TODO: add lbs_NLP, ubs_NLP mutlipliers
    }
    if (get_feasibility_multipliers)
    {
        mem->norm_feas_qp_lam_unslacked_bounds = norm_lam_hard_constr;
        mem->norm_feas_qp_lam_slacked_constraints = norm_lam_slacked_constraints;
    }
    else
    {
        mem->norm_opt_qp_lam_unslacked_bounds = norm_lam_hard_constr;
        mem->norm_opt_qp_lam_slacked_constraints = norm_lam_slacked_constraints;
    }
    // assert(norm_lam_slacked_constraints <= 1.0 + 1e-8); // Slacked multipliers should be in [0,1]
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



/*
Calculates the QP l1 infeasibility based on the additional slack variable values in QP
*/
static double calculate_slacked_qp_l1_infeasibility_from_slacks(ocp_nlp_dims *dims, ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_sqp_wfqp_opts* opts, ocp_qp_out *qp_out)
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
            // Add lower slack
            tmp1 = BLASFEO_DVECEL(qp_out->ux + i, nx[i]+nu[i]+ns[i] + j);
            l1_inf += fmax(0.0, tmp1);

            // Add upper slack
            tmp2 = BLASFEO_DVECEL(qp_out->ux + i, nx[i]+nu[i]+2*ns[i]+nns[i] + j);
            l1_inf += fmax(0.0, tmp2);
            // int index = mem->idxns[i][j];
            // printf("slacks for constraint %d %d: l %e\t u %e\n", i, j, tmp1, tmp2);
        }
    }
    assert(l1_inf > -opts->tol_ineq);
    return l1_inf;
}



/*
This function calculates the l1 infeasibility by calculating the matrix vector product of the
constraints
*/
static double calculate_slacked_qp_l1_infeasibility_manually(ocp_nlp_dims *dims, ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_sqp_wfqp_workspace *work, ocp_nlp_sqp_wfqp_opts* opts, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int *nns = mem->nns;

    int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng;

    double l1_inf = 0.0;
    int i, j;
    double tmp, tmp_bound, mask_value;

    for (i = 0; i <= N; i++)
    {
        // bounds on states and controls
        for (j=0; j<nb[i]; ++j)
        {
            tmp = BLASFEO_DVECEL(qp_out->ux+i, qp_in->idxb[i][j]);
            blasfeo_dvecin1(tmp, &work->nlp_work->tmp_2ni, j);
            blasfeo_dvecin1(tmp, &work->nlp_work->tmp_2ni, nb[i]+ng[i]+j);
        }
        // general linear / linearized!
        // tmp_ni = D * u + C * x
        // lower bounds --> this seems to be correct and in accordance with slack variables
        blasfeo_dgemv_t(nu[i]+nx[i], ng[i], 1.0, qp_in->DCt+i, 0, 0, qp_out->ux+i, 0,
                        0.0, qp_in->d+i, nb[i], &work->nlp_work->tmp_2ni, nb[i]);
        blasfeo_dveccp(ng[i], &work->nlp_work->tmp_2ni, nb[i], &work->nlp_work->tmp_2ni, 2*nb[i]+ng[i]);

        // add slack contributions
        // d[nb:nb+ng] += slack[idx]
        // qp_in->idxs_rev
        if (ns[i]>0)
        {
            for (j = 0; j < nb[i]+ng[i]; j++)
            {
                int slack_index = qp_in->idxs_rev[i][j];
                // maybe we need <=?
                if (slack_index >= 0 && slack_index < ns[i])
                {
                    // add slack contribution for lower and upper constraint
                    // lower
                    BLASFEO_DVECEL(&work->nlp_work->tmp_2ni, j) +=
                            BLASFEO_DVECEL(qp_out->ux+i, slack_index+nx[i]+nu[i]);
                    // upper
                    BLASFEO_DVECEL(&work->nlp_work->tmp_2ni, j+nb[i]+ng[i]) -=
                            BLASFEO_DVECEL(qp_out->ux+i, slack_index+nx[i]+nu[i]+ns[i]+nns[i]);
                }
            }
        }

        // upper bounds (seems to be correct but I do not understand why??)
        // the sign of upper bound d is wrong!! We should use -d. Why is that?
        // blasfeo_dgemv_t(nu[i]+nx[i], ng[i], 1.0, qp_in->DCt+i, 0, 0, qp_out->ux+i, 0,
        //                 0.0, qp_in->d+i, 2*nb[i]+ng[i], &work->nlp_work->tmp_2ni, 2*nb[i]+ng[i]);
        for (j=0; j<2*nb[i]+2*ng[i]; ++j)
        {
            mask_value = BLASFEO_DVECEL(qp_in->d_mask+i, j);
            if (mask_value == 1.0)
            {
                tmp = BLASFEO_DVECEL(&work->nlp_work->tmp_2ni, j);
                tmp_bound = BLASFEO_DVECEL(qp_in->d+i, j);
                if (j < nb[i] + ng[i])
                {
                    // maximum(0, lower_bound - value)
                    l1_inf += fmax(0.0, tmp_bound-tmp);
                    // printf("lower bounds: bound: %.4e, value: %.4e, result: %.4e\n", tmp_bound, tmp, fmax(0.0, tmp_bound-tmp));
                }
                else
                {
                    // upper bounds have the wrong sign!
                    // it is lower_bounds <= value <= -upper_bounds, therefore plus below
                    // printf("upper bounds: value: %.4e, value: %.4e, result: %.4e\n", tmp_bound, tmp, fmax(0.0, tmp_bound+tmp));
                    l1_inf += fmax(0.0, tmp_bound+tmp);
                }
            }
        }
    }
    assert(l1_inf > -opts->tol_ineq);
    return l1_inf;
}



static double calculate_slacked_qp_l1_infeasibility(ocp_nlp_dims *dims, ocp_nlp_sqp_wfqp_memory *mem,
                                                    ocp_nlp_sqp_wfqp_workspace *work, ocp_nlp_sqp_wfqp_opts* opts,
                                                    ocp_qp_in *qp_in, ocp_qp_out *qp_out,
                                                    bool use_slacks)
{
    double l1_inf = 0.0;
    // this is only possible if directly after a QP was solved
    if (use_slacks)
    {
        // seems to be inaccurate. Results are worse!
        l1_inf = calculate_slacked_qp_l1_infeasibility_from_slacks(dims, mem, opts, qp_out);
    }
    else
    {
        l1_inf = calculate_slacked_qp_l1_infeasibility_manually(dims, mem, work, opts, qp_in, qp_out);
    }
    return l1_inf;
}



static void set_non_slacked_l1_penalties(ocp_nlp_config *config, ocp_nlp_dims *dims,
    ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_sqp_wfqp_memory *mem,
    ocp_nlp_workspace *work)
{
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int *nns = mem->nns;
    ocp_qp_in *relaxed_qp_in = mem->relaxed_qp_in;

    // be aware of rqz_QP = [r, q, zl_NLP, zl_QP, zu_NLP, zu_QP]
    for (int stage = 0; stage <= dims->N; stage++)
    {
        // zl_QP
        blasfeo_dvecse(nns[stage], 1.0, relaxed_qp_in->rqz+stage, nu[stage]+nx[stage]+ns[stage]);
        // zu_QP
        blasfeo_dvecse(nns[stage], 1.0, relaxed_qp_in->rqz+stage, nu[stage]+nx[stage]+2*ns[stage]+nns[stage]);
    }
}



static void set_non_slacked_l2_penalties(ocp_nlp_config *config, ocp_nlp_dims *dims,
    ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_sqp_wfqp_memory *mem,
    ocp_nlp_workspace *work)
{
    int *ns = dims->ns;
    int *nns = mem->nns;
    ocp_qp_in *qp_in = mem->nlp_mem->qp_in;
    ocp_qp_in *relaxed_qp_in = mem->relaxed_qp_in;

    // be aware of rqz_QP = [r, q, zl_NLP, zl_QP, zu_NLP, zu_QP]
    for (int stage = 0; stage <= dims->N; stage++)
    {
        // // zu_NLP shift back
        blasfeo_dveccp(ns[stage], qp_in->Z+stage, ns[stage], relaxed_qp_in->Z+stage, ns[stage]+nns[stage]);
        blasfeo_dveccp(ns[stage], qp_in->Z+stage, ns[stage], relaxed_qp_in->Z+stage, ns[stage]+nns[stage]); //-->this was changed!

        // // TODO: rethink! I think we want all l2 slack contributitions to be zero!!
        // // zl_QP
        // // blasfeo_dvecse(nns[stage], 0.0, qp_in->Z+stage, ns[stage]);
        blasfeo_dvecse(nns[stage], 0.0, relaxed_qp_in->Z+stage, ns[stage]);
        // // zu_QP
        // // blasfeo_dvecse(nns[stage], 0.0, qp_in->Z+stage, 2*ns[stage]+nns[stage]);
        blasfeo_dvecse(nns[stage], 0.0, relaxed_qp_in->Z+stage, 2*ns[stage]+nns[stage]);

        // set all to zero?
        // blasfeo_dvecse(2*ns[stage]+2*nns[stage], 0.0, relaxed_qp_in->Z+stage, 0);

        // printf("qp_in->Z %d\n", stage);
        // blasfeo_print_exp_tran_dvec(2*(nns[stage]+ns[stage]), qp_in->Z+stage, 0);
    }
}



static void set_feasibility_multipliers(ocp_nlp_dims *dims,
                            ocp_nlp_sqp_wfqp_memory *mem,
                            ocp_nlp_out *nlp_out)
{
    int *ni = dims->ni;
    int *nx = dims->nx;
    int N = dims->N;

    for (int i=0; i<dims->N; ++i)
    {
        blasfeo_dveccp(2*ni[i], nlp_out->lam+i, 0, mem->lam_feasibility+i, 0);
        if (i < N)
        {
            blasfeo_dveccp(nx[i+1], nlp_out->pi+i, 0, mem->pi_feasibility+i, 0);
        }
    }
}



static double compute_gradient_directional_derivative(ocp_nlp_sqp_wfqp_memory* mem, ocp_nlp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    // Compute the directional derivative of the user-specified objective in the direction qp_out
    double dir_der = 0.0;
    int i, nux;
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int *nns = mem->nns;
    // Sum over stages 0 to N
    for (i = 0; i <= N; i++)
    {
        nux = nx[i] + nu[i];
        // Calculate g.T d
        dir_der += blasfeo_ddot(nux, &qp_out->ux[i], 0, &qp_in->rqz[i], 0);

        // Calculate gradient of slacks.T d_slacks
        // First part of slacks
        dir_der += blasfeo_ddot(ns[i], &qp_out->ux[i], nux, &qp_in->rqz[i], nux);
        // Second part of slacks
        dir_der += blasfeo_ddot(ns[i], &qp_out->ux[i], nux+ns[i]+nns[i], &qp_in->rqz[i], nux+ns[i]+nns[i]);
    }
    return dir_der;
}



static double compute_qp_objective_value(ocp_nlp_sqp_wfqp_memory* mem, ocp_nlp_dims *dims,
                                  ocp_qp_in *qp_in,
                                  ocp_qp_out *qp_out,
                                  ocp_nlp_workspace *nlp_work)
{
    // Compute the QP objective function value corresponding to the user specified objective.
    double qp_cost = 0.0;
    int i, nux, ns, nns;
    int N = dims->N;
    // Sum over stages 0 to N
    for (i = 0; i <= N; i++)
    {
        nux = dims->nx[i] + dims->nu[i];
        ns = dims->ns[i];
        nns = mem->nns[i];
        // Calculate 0.5 * d.T H d
        blasfeo_dsymv_l(nux, 0.5, &qp_in->RSQrq[i], 0, 0, &qp_out->ux[i], 0,
                        0.0, &qp_out->ux[i], 0, &nlp_work->tmp_nv, 0);
        qp_cost += blasfeo_ddot(nux, &qp_out->ux[i], 0, &nlp_work->tmp_nv, 0);

        // slack QP objective value, compare to computation in cost modules;
        // lower
        // tmp_nv = 2 * z + Z .* slack;
        blasfeo_dveccpsc(ns, 2.0, &qp_out->ux[i], nux, &nlp_work->tmp_nv, 0);
        blasfeo_dvecmulacc(ns, &qp_in->Z[i], 0, &qp_out->ux[i], nux, &nlp_work->tmp_nv, 0);
        // qp_cost += .5 * (tmp_nv .* slack)
        qp_cost += 0.5 * blasfeo_ddot(ns, &nlp_work->tmp_nv, 0, &qp_out->ux[i], nux);

        // upper
        // tmp_nv = 2 * z + Z .* slack;
        blasfeo_dveccpsc(ns, 2.0, &qp_out->ux[i], nux+ns+nns, &nlp_work->tmp_nv, 0);
        blasfeo_dvecmulacc(ns, &qp_in->Z[i], 0, &qp_out->ux[i], nux+ns+nns, &nlp_work->tmp_nv, 0);
        // qp_cost += .5 * (tmp_nv .* slack)
        qp_cost += 0.5 * blasfeo_ddot(ns, &nlp_work->tmp_nv, 0, &qp_out->ux[i], nux+ns+nns);

        // Calculate g.T d
        qp_cost += blasfeo_ddot(nux, &qp_out->ux[i], 0, &qp_in->rqz[i], 0);

        // Calculate gradient of slacks.T d_slacks
        // TODO: either comment or formula is not correct, this computs slack.T * z
        // lower of slacks
        qp_cost += blasfeo_ddot(ns, &qp_out->ux[i], nux, &qp_in->rqz[i], nux);
        // upper of slacks
        qp_cost += blasfeo_ddot(ns, &qp_out->ux[i], nux+ns+nns, &qp_in->rqz[i], nux+ns+nns);
    }
    return qp_cost;
}



static void print_indices(ocp_nlp_dims *dims, ocp_nlp_sqp_wfqp_workspace *work, ocp_nlp_sqp_wfqp_memory *mem)
{
    ocp_nlp_workspace *nlp_work = work->nlp_work;
    int ni, ns, nns;
    int *idxs = nlp_work->tmp_nins;
    for (int stage = 0; stage <= dims->N; stage++)
    {
        ns = dims->ns[stage];
        ni = dims->ni[stage];
        nns = mem->nns[stage];
        printf("stage %d: ni %d ns %d nns %d\n", stage, ni, ns, nns);
        printf("got idxs at stage %d\n", stage);
        for (int i=0; i<ns; i++)
            printf("%d ", idxs[i]);
        printf("\n");

        printf("got idxns at stage %d\n", stage);
        for (int i=0; i<nns; i++)
            printf("%d ", mem->idxns[stage][i]);
        printf("\n");
    }
}



/************************************************
 * functions for QP preparation
 ************************************************/
// TODO: we still need this somewhere such that we can use the feasibility QP with constraint Hessian matrix
// static void ocp_nlp_sqp_wfqp_prepare_hessian_evaluation(ocp_nlp_config *config,
//     ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
//     ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_workspace *work)
// {
//     ocp_nlp_memory *nlp_mem = mem->nlp_mem;

//     int N = dims->N;

//     // int *nv = dims->nv;
//     int *nx = dims->nx;
//     int *nu = dims->nu;

//     if (!nlp_mem->compute_hess)
//     {
//         printf("ocp_nlp_sqp_wfqp_prepare_hessian_evaluation: constant hessian not supported!\n\n");
//         exit(1);
//     }

// #if defined(ACADOS_WITH_OPENMP)
//     #pragma omp parallel for
// #endif
//     for (int i = 0; i <= N; i++)
//     {
//         // TODO: first compute cost hessian (without adding) and avoid setting everything to zero?
//         // init Hessians to 0

//         // TODO: avoid setting qp_in->RSQ to zero in ocp_nlp_approximate_qp_matrices?
//         blasfeo_dgese(nu[i] + nx[i], nu[i] + nx[i], 0.0, mem->RSQ_constr+i, 0, 0);
//         blasfeo_dgese(nu[i] + nx[i], nu[i] + nx[i], 0.0, mem->RSQ_cost+i, 0, 0);
//     }

//     for (int i = 0; i < N; i++)
//     {
//         config->dynamics[i]->memory_set_RSQrq_ptr(mem->RSQ_constr+i, nlp_mem->dynamics[i]);
//     }
//     for (int i = 0; i <= N; i++)
//     {
//         config->cost[i]->memory_set_RSQrq_ptr(mem->RSQ_cost+i, nlp_mem->cost[i]);
//         config->constraints[i]->memory_set_RSQrq_ptr(mem->RSQ_constr+i, nlp_mem->constraints[i]);
//     }
//     return;
// }

// update QP rhs for SQP (step prim var, abs dual var)
// - use cost gradient and dynamics residual from memory
// - evaluate constraints wrt bounds -> allows to update all bounds between preparation and feedback phase.
void ocp_nlp_sqp_wfqp_approximate_qp_constraint_vectors(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
    ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_workspace *work, int sqp_iter)
{
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_qp_in *relaxed_qp_in = mem->relaxed_qp_in;
    int N = dims->N;

    int *ns = dims->ns;
    int *ni = dims->ni;
    int *nns = mem->nns;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // b is set by pointer of nominal QP

        // evaluate constraint residuals
        config->constraints[i]->update_qp_vectors(config->constraints[i], dims->constraints[i],
            in->constraints[i], opts->constraints[i], nlp_mem->constraints[i], work->constraints[i]);

        // copy ineq function value into nlp mem, then into QP
        struct blasfeo_dvec *ineq_fun = config->constraints[i]->memory_get_fun_ptr(nlp_mem->constraints[i]);
        blasfeo_dveccp(2 * ni[i], ineq_fun, 0, nlp_mem->ineq_fun + i, 0);

        // d
        int n_nominal_ineq_nlp = ni[i] - ns[i];

        // blasfeo_dveccp(2 * ni[i], nlp_mem->ineq_fun + i, 0, relaxed_qp_in->d + i, 0);
        blasfeo_dveccp(2*n_nominal_ineq_nlp+ns[i], nlp_mem->ineq_fun + i, 0, relaxed_qp_in->d + i, 0);
        blasfeo_dveccp(ns[i], nlp_mem->ineq_fun + i, 2*n_nominal_ineq_nlp+ns[i], relaxed_qp_in->d + i, 2*n_nominal_ineq_nlp+ns[i]+nns[i]);
    }
    // setup d_mask
    if (sqp_iter == 0)
    {
        int offset_dmask;
        for (int i=0; i<=dims->N; i++)
        {
            offset_dmask = 2*(dims->nb[i]+dims->ng[i]+dims->ni_nl[i]);
            blasfeo_dveccp(offset_dmask, nlp_mem->qp_in->d_mask+i, 0, relaxed_qp_in->d_mask+i, 0);
            blasfeo_dvecse(2*relaxed_qp_in->dim->ns[i], 1.0, relaxed_qp_in->d_mask+i,offset_dmask);
        }
    }
}



static void ocp_nlp_sqp_wfqp_setup_feasibility_qp_objective(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_sqp_wfqp_opts *opts,
    ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_workspace *work)
{
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    // ocp_qp_in *qp_in = nlp_mem->qp_in;
    ocp_qp_in *qp_in = mem->relaxed_qp_in;
    // TODO: ocp_qp_in *qp_in = mem->relaxed_qp_in;
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
        // NOTE: I think this TODO is from before we switched to using objective multiplier and it should be fine now.

        // if (opts->use_exact_hessian_in_feas_qp)
        // {
        //     // Either we use the exact objective Hessian
        //     blasfeo_dgecp(nxu, nxu, mem->RSQ_constr+i, 0, 0, qp_in->RSQrq+i, 0, 0);
        //     blasfeo_dgead(nxu, nxu, objective_multiplier, mem->RSQ_cost+i, 0, 0, qp_in->RSQrq+i, 0, 0);
        // }
        // else
        // {
        // We use the identity matrix Hessian
        blasfeo_dgese(nxu, nxu, 0.0, qp_in->RSQrq+i, 0, 0);
        blasfeo_ddiare(nxu, 1e-4, qp_in->RSQrq+i, 0, 0);
        // indicate that Hessian is diagonal
        qp_in->diag_H_flag[i] = 1;
        // }

        // Z -- slack matrix
        blasfeo_dveccpsc(ns[i], 0.0, mem->Z_cost_module+i, 0, qp_in->Z+i, 0);
        blasfeo_dveccpsc(ns[i], 0.0, mem->Z_cost_module+i, 0, qp_in->Z+i, ns[i]+mem->nns[i]);

        /* vectors */
        // g
        blasfeo_dveccpsc(nx[i]+nu[i]+ns[i], 0.0, nlp_mem->cost_grad + i, 0, qp_in->rqz + i, 0);
        blasfeo_dveccpsc(ns[i], 0.0, nlp_mem->cost_grad + i, nx[i]+nu[i]+ns[i], qp_in->rqz + i, nx[i]+nu[i]+ns[i]+nns[i]);

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
    ocp_nlp_timings *nlp_timings = nlp_mem->nlp_timings;
    int qp_status = ACADOS_SUCCESS;

    // (typically) no warm start at first iteration
    if (sqp_iter == 0 && !opts->warm_start_first_qp)
    {
        int tmp_int = 0;
        qp_solver->opts_set(qp_solver, nlp_opts->qp_solver_opts, "warm_start", &tmp_int);
    }

    // ocp_nlp_add_levenberg_marquardt_term(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, mem->alpha, sqp_iter);
    // regularize Hessian
    acados_tic(&timer1);
    config->regularize->regularize(config->regularize, dims->regularize, nlp_opts->regularize, nlp_mem->regularize_mem);
    nlp_timings->time_reg += acados_toc(&timer1);
    // Show input to QP
    if (nlp_opts->print_level > 3)
    {
        printf("\n\nSQP: ocp_qp_in at iteration %d\n", sqp_iter);
        print_ocp_qp_dims(qp_in->dim);
        print_ocp_qp_in(qp_in);
    }

#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
    ocp_nlp_dump_qp_in_to_file(qp_in, sqp_iter, 0);
#endif

    if (solve_feasibility_qp)
    {
        qp_status = ocp_nlp_solve_qp_and_correct_dual(config, dims, nlp_opts,
                                                    nlp_mem, nlp_work, false,
                                                    qp_in, qp_out, &mem->relaxed_qp_solver);
    }
    else
    {
        qp_status = ocp_nlp_solve_qp_and_correct_dual(config, dims, nlp_opts,
                                                    nlp_mem, nlp_work, false,
                                                    NULL, NULL, NULL);
    }

    // restore default warm start
    if (sqp_iter==0)
    {
        qp_solver->opts_set(qp_solver, nlp_opts->qp_solver_opts, "warm_start", &opts->qp_warm_start);
    }

    if (nlp_opts->print_level > 3)
    {
        printf("\n\nSQP: ocp_qp_out at iteration %d\n", sqp_iter);
        print_ocp_qp_dims(qp_out->dim);
        print_ocp_qp_out(qp_out);
    }

#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
    ocp_nlp_dump_qp_out_to_file(qp_out, sqp_iter, 0);
#endif

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
    return qp_status;
}

/************************************************
* Byrd-Omojokun Subproblem Stuff
************************************************/
/*
Adjusts the bounds of the QP with the given slack variables, such that the
resulting QP has always a feasible solution.
*/
static void setup_byrd_omojokun_bounds(ocp_nlp_dims *dims, ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_sqp_wfqp_workspace *work, ocp_nlp_sqp_wfqp_opts* opts, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int *nns = mem->nns;

    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ni_nl = dims->ni_nl;

    int i, j;
    double tmp_lower, tmp_upper;

    for (i = 0; i <= N; i++)
    {
        // printf("i = %d\n", i);
        // printf("ns[i] = %d\n", dims->ns[i]);
        // printf("nns[i] = %d\n", mem->nns[i]);
        for (j=0; j<mem->nns[i]; ++j)
        {
            // printf("j = %d\n", j);
            int constr_index = mem->idxns[i][j]; //
            int slack_index = mem->relaxed_qp_in->idxs_rev[i][constr_index];
            // printf("slack_index = %d\n", slack_index);
            // get lower slack
            tmp_lower = BLASFEO_DVECEL(qp_out->ux + i, nx[i]+nu[i]+slack_index);
            // lower_bound - value
            // TODO: d should be copied such that we can extract in python
            BLASFEO_DVECEL(qp_in->d+i, constr_index) -= tmp_lower;

            // get upper slack
            tmp_upper = BLASFEO_DVECEL(qp_out->ux + i, nx[i]+nu[i]+ns[i]+nns[i] + slack_index);
            // upper bounds have the wrong sign!
            // it is lower_bounds <= value <= -upper_bounds, therefore plus below
            // for the slacks with upper bound we have value - slack, therefore
            // value <= -upper_bound + slack,
            // therefore we store upper_bound - slack??
            BLASFEO_DVECEL(qp_in->d+i, nb[i] + ng[i] + ni_nl[i] + constr_index) -= tmp_upper;
        }
    }
}

/*
Solve feasibility QP first and then solve standard QP with adjusted bounds
*/
static int byrd_omojokun_direction_computation(ocp_nlp_dims *dims,
                                            ocp_nlp_config *config,
                                            ocp_nlp_sqp_wfqp_opts *opts,
                                            ocp_nlp_opts *nlp_opts,
                                            ocp_nlp_in *nlp_in,
                                            ocp_nlp_out *nlp_out,
                                            ocp_nlp_sqp_wfqp_memory *mem,
                                            ocp_nlp_sqp_wfqp_workspace *work,
                                            double current_l1_infeasibility,
                                            int sqp_iter,
                                            acados_timer timer0,
                                            acados_timer timer1)
{
    ocp_nlp_memory* nlp_mem = mem->nlp_mem;
    ocp_nlp_workspace* nlp_work = work->nlp_work;
    ocp_qp_in *nominal_qp_in = nlp_mem->qp_in;
    ocp_qp_out *nominal_qp_out = nlp_mem->qp_out;

    ocp_qp_in *relaxed_qp_in = mem->relaxed_qp_in;
    ocp_qp_out *relaxed_qp_out = mem->relaxed_qp_out;

    ocp_nlp_timings *nlp_timings = nlp_mem->nlp_timings;
    qp_info* qp_info_;

    int qp_status;
    int qp_iter = 0;
    double pred_l1_inf_QP_feasibility;
    double l1_inf_QP_feasibility;

    /* Solve Feasibility QP: Only gradient of slack variables */
    print_debug_output("Solve Feasibility QP!\n", nlp_opts->print_level, 2);
    // print_ocp_qp_dims(relaxed_qp_in->dim);
    // print_ocp_qp_in(relaxed_qp_in);
    qp_status = prepare_and_solve_QP(config, opts, relaxed_qp_in, relaxed_qp_out, dims, mem, nlp_in, nlp_out,
                nlp_mem, nlp_work, sqp_iter, true, timer0, timer1);
    ocp_qp_out_get(relaxed_qp_out, "qp_info", &qp_info_);
    qp_iter += qp_info_->num_iter;
    if (qp_status != ACADOS_SUCCESS)
    {
        if (nlp_opts->print_level >=1)
        {
            printf("\n Error in Feasibility QP in iteration %d, got qp_status %d!\n", qp_iter, qp_status);
        }
        nlp_mem->status = ACADOS_QP_FAILURE;
        nlp_mem->iter = sqp_iter;
        nlp_timings->time_tot = acados_toc(&timer0);
        // #if defined(ACADOS_WITH_OPENMP)
        //         // restore number of threads
        //         int num_threads_bkp = omp_get_num_threads();
        //         omp_set_num_threads(num_threads_bkp);
        // #endif
        return nlp_mem->status;
    }
    compute_qp_multiplier_norm_inf(mem, dims, relaxed_qp_out, relaxed_qp_in, true);

    l1_inf_QP_feasibility = calculate_slacked_qp_l1_infeasibility(dims, mem, work, opts, relaxed_qp_in, relaxed_qp_out, opts->use_QP_l1_inf_from_slacks);
    pred_l1_inf_QP_feasibility = calculate_predicted_l1_inf_reduction(opts, current_l1_infeasibility, l1_inf_QP_feasibility);
    mem->pred_l1_inf_QP_optimality = pred_l1_inf_QP_feasibility;

    print_debug_output_double("Feas QP: Multiplier norm: pi", mem->norm_feas_qp_pi, nlp_opts->print_level, 2);
    print_debug_output_double("Feas QP: Multiplier norm: lam slacked", mem->norm_feas_qp_lam_slacked_constraints, nlp_opts->print_level, 2);
    print_debug_output_double("Feas QP: Multiplier norm: lam unslacked", mem->norm_opt_qp_lam_unslacked_bounds, nlp_opts->print_level, 2);
    print_debug_output_double("Feas QP: l1_inf_feas: ", l1_inf_QP_feasibility, nlp_opts->print_level, 2);
    print_debug_output_double("Feas QP: pred_l1_inf_QP: ", pred_l1_inf_QP_feasibility, nlp_opts->print_level, 2);
    // assert(pred_l1_inf_QP_feasibility > -1e2*opts->tol_ineq);

    /* Solve the nominal QP with updated bounds*/
    print_debug_output("Solve Nominal QP!\n", nlp_opts->print_level, 2);
    setup_byrd_omojokun_bounds(dims, mem, work, opts, nominal_qp_in, relaxed_qp_out);
    // solve_feasibility_qp = false below
    qp_status = prepare_and_solve_QP(config, opts, nominal_qp_in, nominal_qp_out, dims, mem, nlp_in, nlp_out,
                                     nlp_mem, nlp_work, sqp_iter, false, timer0, timer1);
    if (qp_status != ACADOS_SUCCESS)
    {
        if (nlp_opts->print_level >=1)
        {
            printf("\n Error in Nominal QP in iteration %d, got qp_status %d!\n", qp_iter, qp_status);
        }
        nlp_mem->status = ACADOS_QP_FAILURE;
        nlp_mem->iter = sqp_iter;
        nlp_timings->time_tot = acados_toc(&timer0);
        // #if defined(ACADOS_WITH_OPENMP)
        //         // restore number of threads
        //         int num_threads_bkp = omp_get_num_threads();
        //         omp_set_num_threads(num_threads_bkp);
        // #endif
        return nlp_mem->status;
    }
    compute_qp_multiplier_norm_inf(mem, dims, nominal_qp_out, nominal_qp_in, false);
    return qp_status;
}

/*
Depending on the search direction mode, a search direction is calculated
*/
static int calculate_search_direction(ocp_nlp_dims *dims,
                                        ocp_nlp_config *config,
                                        ocp_nlp_sqp_wfqp_opts *opts,
                                        ocp_nlp_opts *nlp_opts,
                                        ocp_nlp_in *nlp_in,
                                        ocp_nlp_out *nlp_out,
                                        ocp_nlp_sqp_wfqp_memory *mem,
                                        ocp_nlp_sqp_wfqp_workspace *work,
                                        double current_l1_infeasibility,
                                        int sqp_iter,
                                        acados_timer timer0,
                                        acados_timer timer1)
{
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    qp_info* qp_info_;
    int qp_iter = 0;
    int search_direction_status;
    int solved_nominal_before = 0;
    if (mem->search_direction_mode == NOMINAL_QP)
    {
        // if the QP can be solved and the status is good, we
        // return 0
        // otherwise, we change the mode to Byrd-Omojokun and we continue
        // below!
        search_direction_status = prepare_and_solve_QP(config, opts, nlp_mem->qp_in, nlp_mem->qp_out,
                                                       dims, mem, nlp_in, nlp_out, nlp_mem, work->nlp_work,
                                                       sqp_iter, false, timer0, timer1);
                                                       ocp_qp_out_get(nlp_mem->qp_out, "qp_info", &qp_info_);
        qp_iter += qp_info_->num_iter;
        if (search_direction_status != ACADOS_SUCCESS)
        {
            if (nlp_opts->print_level >=1)
            {
                printf("\n Error in Nominal QP in iteration %d, got qp_status %d!\n", qp_iter, search_direction_status);
                printf("Switch to Byrd-Omojokun mode!\n");
            }
            mem->search_direction_mode = BYRD_OMOJOKUN;
            solved_nominal_before = 1;
        }
        else
        {
            compute_qp_multiplier_norm_inf(mem, dims, nlp_mem->qp_out, nlp_mem->qp_in, false);
            mem->search_direction_type = "N";
            return ACADOS_SUCCESS;
        }
    }
    if (mem->search_direction_mode == BYRD_OMOJOKUN)
    {
        // We solve two QPs and return the search direction that we found!
        // if the second QP is feasible, we change back to nominal QP mode.
        // Maybe we want some kind of watchdog, if for two consecutive QPs this holds
        // then we switch back
        search_direction_status = byrd_omojokun_direction_computation(dims,
                                                                    config,
                                                                    opts,
                                                                    nlp_opts,
                                                                    nlp_in,
                                                                    nlp_out,
                                                                    mem,
                                                                    work,
                                                                    current_l1_infeasibility,
                                                                    sqp_iter,
                                                                    timer0,
                                                                    timer1);
        // TODO: solve this below!!!
        // search_direction_status = 0;
        if (solved_nominal_before)
        {
            mem->search_direction_type = "NFN";
        }
        else
        {
            mem->search_direction_type = "FN";
        }
        // Here is something missing that we might switch back to the nominal mode
        double l1_inf = calculate_slacked_qp_l1_infeasibility_from_slacks(dims, mem, opts, mem->relaxed_qp_out);
        printf("Slack sum: %10.4e\n", l1_inf);
        if (l1_inf/(fmax(1.0, (double) mem->absolute_nns)) < opts->tol_ineq)
        {
            mem->watchdog_zero_slacks_counter += 1;
        }

        if (opts->allow_direction_mode_switch && mem->watchdog_zero_slacks_counter == opts->watchdog_zero_slacks_max)
        {
            mem->watchdog_zero_slacks_counter = 0;
            mem->search_direction_mode = NOMINAL_QP;
        }

        return search_direction_status;
        // Maybe we want to switch to a full feasibility restoration phase
        // if the NLP seems to be infeasible?
    }
    else if (mem->search_direction_mode == FEASIBILITY_QP)
    {
        // for the moment we do nothing here!
        printf("Feasibility mode not implemented at the moment!\n");
        mem->search_direction_type = "F";
        return 1;
    }
    else
    {
        printf("Wrong search direction mode\n");
        return 1;
    }
}



/********************************
* Functions for relaxed QP
*********************************/
static void set_relaxed_qp_in_matrix_pointers(ocp_nlp_sqp_wfqp_memory *mem, ocp_qp_in *qp_in)
{
    mem->relaxed_qp_in->BAbt = qp_in->BAbt; // dynamics matrix & vector work space
    mem->relaxed_qp_in->b = qp_in->b; // dynamics vector work space
    // mem->relaxed_qp_in->RSQrq -->  Hessians should be different Hessian
    mem->relaxed_qp_in->DCt = qp_in->DCt; // inequality constraints matrix
    mem->relaxed_qp_in->idxb = qp_in->idxb;
    mem->relaxed_qp_in->idxe = qp_in->idxe;
    // mem->relaxed_qp_in->d_mask --> different due to slacks
    // mem->relaxed_qp_in->diag_H_flag --> if identity Hessian used in feasibility QP, flag set elsewhere
    // mem->relaxed_qp_in->m = qp_in->m; // TODO: Not sure what happens here
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
    ocp_nlp_res *nlp_res = nlp_mem->nlp_res;
    ocp_nlp_timings *nlp_timings = nlp_mem->nlp_timings;

    ocp_nlp_sqp_wfqp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    ocp_qp_in *nominal_qp_in = nlp_mem->qp_in;
    ocp_qp_out *nominal_qp_out = nlp_mem->qp_out;
    ocp_qp_in *relaxed_qp_in = mem->relaxed_qp_in;
    ocp_qp_out *relaxed_qp_out = mem->relaxed_qp_out;

    // zero timers
    ocp_nlp_timings_reset(nlp_timings);

    int qp_status = 0;
    int qp_iter = 0;
    mem->alpha = 0.0;
    mem->step_norm = 0.0;
    mem->nlp_mem->status = ACADOS_READY;
    mem->search_direction_type = "-";
    nlp_mem->objective_multiplier = 1.0;
    mem->search_direction_mode = opts->search_direction_mode;
    mem->watchdog_zero_slacks_counter = 0;

#if defined(ACADOS_WITH_OPENMP)
    // backup number of threads
    int num_threads_bkp = omp_get_num_threads();
    // set number of threads
    omp_set_num_threads(opts->nlp_opts->num_threads);
#endif

    mem->absolute_nns = 0;
    for (int i = 0; i <= dims->N; i++)
    {
        mem->absolute_nns += mem->nns[i];
    }

    set_relaxed_qp_in_matrix_pointers(mem, nominal_qp_in);

    ocp_nlp_initialize_submodules(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
    // feasibility QP objective always constant! Only done once!
    ocp_nlp_sqp_wfqp_setup_feasibility_qp_objective(config, dims, nlp_in, nlp_out, opts, mem, nlp_work);

    set_non_slacked_l2_penalties(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work);
    set_non_slacked_l1_penalties(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work);
    set_feasibility_multipliers(dims, mem, nlp_out);

    /************************************************
     * main sqp loop
     ************************************************/
    int sqp_iter = 0;
    double prev_levenberg_marquardt = 0.0;
    int search_direction_status = 0;
    double pred_l1_inf_search_direction;
    double l1_inf_search_direction;

    if (nlp_opts->print_level > 1)
    {
        print_indices(dims, work, mem);
    }

    for (; sqp_iter <= opts->max_iter; sqp_iter++) // <= needed such that after last iteration KKT residuals are checked before max_iter is thrown.
    {
        // We always evaluate the residuals until the last iteration
        // If the option "eval_residual_at_max_iter" is set, we also
        // evaluate the residuals after the last iteration.
        if (sqp_iter != opts->max_iter || opts->eval_residual_at_max_iter)
        {
            /* Prepare the QP data */
            // linearize NLP and update QP matrices
            // for nominal QP only. relaxed QP has identity Hessian
            // TODO: this seems outdated at the moment!
            // ocp_nlp_sqp_wfqp_prepare_hessian_evaluation(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work);
            acados_tic(&timer1);
            // nominal QP solver
            ocp_nlp_approximate_qp_matrices(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
            ocp_nlp_approximate_qp_vectors_sqp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);

            // relaxed QP solver
            // matrices for relaxed QP solver evaluated in nominal QP solver
            ocp_nlp_sqp_wfqp_approximate_qp_constraint_vectors(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work, sqp_iter);

            if (nlp_opts->with_adaptive_levenberg_marquardt || config->globalization->needs_objective_value() == 1)
            {
                ocp_nlp_get_cost_value_from_submodules(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
            }
            //
            nlp_timings->time_lin += acados_toc(&timer1);
            // compute nlp residuals
            ocp_nlp_res_compute(dims, nlp_opts, nlp_in, nlp_out, nlp_res, nlp_mem, nlp_work);
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

        /* Output */
        if (nlp_opts->print_level > 0)
        {
            print_iteration(sqp_iter, config, nlp_res, mem, nlp_opts, prev_levenberg_marquardt, qp_status, qp_iter);
        }
        prev_levenberg_marquardt = nlp_opts->levenberg_marquardt;

        /* Termination */
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

        double current_l1_infeasibility = ocp_nlp_get_l1_infeasibility(config, dims, nlp_mem);
        print_debug_output_double("Current l1 infeasibility: ", current_l1_infeasibility, nlp_opts->print_level, 2);

        /* search direction computation */
        search_direction_status = calculate_search_direction(dims,
                                                            config,
                                                            opts,
                                                            nlp_opts,
                                                            nlp_in,
                                                            nlp_out,
                                                            mem,
                                                            work,
                                                            current_l1_infeasibility,
                                                            sqp_iter,
                                                            timer0,
                                                            timer1);
        if (search_direction_status != ACADOS_SUCCESS)
        {
            return nlp_mem->status;
        }

        // Log the qp stats. At the moment we sum up the number of total QP iterations
        // The solver anyway terminates if a QP was not solved correctly at this point
        if (sqp_iter+1 < mem->stat_m)
        {
            mem->stat[mem->stat_n*(sqp_iter+1)+4] = qp_status;
            mem->stat[mem->stat_n*(sqp_iter+1)+5] = qp_iter;
        }

        // We want to keep this! Since l1 QP should have the same as our ByrdOmojokun-QP
        l1_inf_search_direction = calculate_slacked_qp_l1_infeasibility(dims, mem, work, opts, relaxed_qp_in, relaxed_qp_out, false);
        pred_l1_inf_search_direction = calculate_predicted_l1_inf_reduction(opts, current_l1_infeasibility, l1_inf_search_direction);
        print_debug_output_double("pred_l1_inf_search_direction: ", pred_l1_inf_search_direction, nlp_opts->print_level, 2);

        // Compute the step norm
        if (opts->tol_min_step_norm > 0.0 || nlp_opts->log_primal_step_norm)
        {
            // For the moment we do not care about the artificial slack variables to keep problem
            // feasible
            mem->step_norm = ocp_qp_out_compute_primal_nrm_inf(nominal_qp_out);
            if (nlp_opts->log_primal_step_norm)
                mem->primal_step_norm[sqp_iter] = mem->step_norm;
        }

        /* globalization */
        // Calculate optimal QP objective (needed for globalization)
        if (config->globalization->needs_qp_objective_value() == 1)
        {
            // nlp_mem->qp_cost_value = compute_qp_objective_value(mem, dims, qp_in, qp_out, nlp_work);
            nlp_mem->qp_cost_value = compute_qp_objective_value(mem, dims, nominal_qp_in, nominal_qp_out, nlp_work);
            nlp_mem->predicted_infeasibility_reduction = pred_l1_inf_search_direction;
            nlp_mem->predicted_optimality_reduction = -compute_gradient_directional_derivative(mem, dims, nominal_qp_in, nominal_qp_out);
            print_debug_output_double("pred_opt_search_direction: ", nlp_mem->predicted_optimality_reduction, nlp_opts->print_level, 2);
        }
        // NOTE on timings: currently all within globalization is accounted for within time_glob.
        //   QP solver times could be also attributed there alternatively. Cleanest would be to save them seperately.
        acados_tic(&timer1);

        int globalization_status;
        globalization_status = config->globalization->find_acceptable_iterate(config, dims, nlp_in,
                                                                              nlp_out, nlp_mem, mem,
                                                                              nlp_work, nlp_opts, &mem->alpha);
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

    nlp_mem->workspace_size = ocp_nlp_workspace_calculate_size(config, dims, opts->nlp_opts, nlp_in);

    ocp_nlp_sqp_wfqp_workspace *work = work_;
    ocp_nlp_sqp_wfqp_cast_workspace(config, dims, opts, nlp_in, mem, work);
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
    }

    // set idxs_rev of relaxed QP
    for (int stage = 0; stage <= dims->N; stage++)
    {
        config->constraints[stage]->model_get(config->constraints[stage], dims->constraints[stage], nlp_in->constraints[stage], "idxs", idxs);
        int *relaxed_idxs_rev = mem->relaxed_qp_in->idxs_rev[stage];
        int *idxns = mem->idxns[stage];
        nns = mem->nns[stage];
        ns = dims->ns[stage];
        for (int i=0; i<ns; i++)
        {
            relaxed_idxs_rev[idxs[i]] = i;
        }
        for (int i=0; i<nns; i++)
        {
            relaxed_idxs_rev[idxns[i]] = i+dims->ns[stage];
        }
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
    else if (!strcmp("relaxed_qp_in", field))
    {
        void **value = return_value_;
        *value = mem->relaxed_qp_in;
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
    // TODO: make sure all functions in ocp_nlp_config are defined!
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
    config->step_update = &ocp_nlp_update_variables_sqp;
    config->is_real_time_algorithm = &ocp_nlp_sqp_wfqp_is_real_time_algorithm;

    return;
}
