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
#include "acados/utils/math.h"
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

    // this first !!!
    ocp_nlp_opts_initialize_default(config, dims, nlp_opts);

    // SQP opts
    opts->nlp_opts->max_iter = 20;

    opts->nlp_opts->eval_residual_at_max_iter = true;
    opts->use_QP_l1_inf_from_slacks = false; // if manual calculation used, results seem more accurate and solver performs better!

    opts->use_constraint_hessian_in_feas_qp = false;
    opts->search_direction_mode = NOMINAL_QP;
    opts->watchdog_zero_slacks_max = 2;
    opts->allow_direction_mode_switch_to_nominal = true;
    opts->feasibility_qp_hessian_scalar = 1e-4;
    opts->log_pi_norm_inf = true;
    opts->log_lam_norm_inf = true;

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
    }
    else // nlp opts
    {
        if (!strcmp(field, "use_constraint_hessian_in_feas_qp"))
        {
            bool* use_constraint_hessian_in_feas_qp = (bool *) value;
            opts->use_constraint_hessian_in_feas_qp = *use_constraint_hessian_in_feas_qp;
        }
        else if (!strcmp(field, "search_direction_mode"))
        {
            bool* search_direction_mode = (bool *) value;
            opts->search_direction_mode = *search_direction_mode;
        }
        else if (!strcmp(field, "allow_direction_mode_switch_to_nominal"))
        {
            bool* allow_direction_mode_switch_to_nominal = (bool *) value;
            opts->allow_direction_mode_switch_to_nominal = *allow_direction_mode_switch_to_nominal;
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
    // relaxed qpscaling
    size += ocp_nlp_qpscaling_memory_calculate_size(dims->qpscaling, nlp_opts->qpscaling, dims->relaxed_qp_solver->orig_dims);

    // stat
    int stat_m = opts->nlp_opts->max_iter+1;
    int stat_n = 13;
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

        // Z_cost_module
        size += blasfeo_memsize_dvec(2*dims->ns[stage]);

        // RSQ_cost, RSQ_constr
        size += 2*blasfeo_memsize_dmat(dims->nx[stage]+dims->nu[stage], dims->nx[stage]+dims->nu[stage]);
    }
    // nns
    size += (N+1) * sizeof(int);
    // Z_cost_module
    size += (N + 1) * sizeof(struct blasfeo_dvec);
    // RSQ_cost, RSQ_constr
    size += 2*(N + 1) * sizeof(struct blasfeo_dmat);

    // search_direction_type
    size += MAX_STR_LEN * sizeof(char);

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
    // qpscaling
    mem->relaxed_qpscaling_mem = ocp_nlp_qpscaling_memory_assign(dims->relaxed_qpscaling, opts->nlp_opts->qpscaling, dims->relaxed_qp_solver->orig_dims, c_ptr);
    c_ptr += ocp_nlp_qpscaling_memory_calculate_size(dims->relaxed_qpscaling, opts->nlp_opts->qpscaling, dims->relaxed_qp_solver->orig_dims);

    // nominal
    mem->relaxed_qp_solver.config = config->relaxed_qp_solver;
    mem->relaxed_qp_solver.dims = dims->relaxed_qp_solver;
    mem->relaxed_qp_solver.opts = opts->nlp_opts->qp_solver_opts;
    mem->relaxed_qp_solver.mem = mem->relaxed_qp_solver_mem;
    mem->relaxed_qp_solver.work = mem->relaxed_qp_solver_work;

    set_relaxed_qp_in_matrix_pointers(mem);

    // Z_cost_module
    assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->Z_cost_module, &c_ptr);

    // RSQ_cost
    assign_and_advance_blasfeo_dmat_structs(N + 1, &mem->RSQ_cost, &c_ptr);
    // RSQ_constr
    assign_and_advance_blasfeo_dmat_structs(N + 1, &mem->RSQ_constr, &c_ptr);

    // stat
    mem->stat_m = opts->nlp_opts->max_iter+1;
    mem->stat_n = 13;
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

    // sum of all nns[i]
    mem->absolute_nns = 0;
    for (int i = 0; i <= dims->N; i++)
    {
        mem->absolute_nns += mem->nns[i];
    }

    mem->nlp_mem->status = ACADOS_READY;
    assign_and_advance_char(MAX_STR_LEN, &mem->search_direction_type, &c_ptr);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);
    // matrices first: RSQ
    for (int i = 0; i <= N; ++i)
    {
        assign_and_advance_blasfeo_dmat_mem(dims->nx[i]+dims->nu[i], dims->nx[i]+dims->nu[i], mem->RSQ_cost + i, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(dims->nx[i]+dims->nu[i], dims->nx[i]+dims->nu[i], mem->RSQ_constr + i, &c_ptr);
    }
    // Z_cost_module
    for (int i = 0; i <= N; ++i)
    {
        assign_and_advance_blasfeo_dvec_mem(2*dims->ns[i], mem->Z_cost_module + i, &c_ptr);
    }
    assert((char *) raw_memory + ocp_nlp_sqp_wfqp_memory_calculate_size(config, dims, opts, in) >= c_ptr);

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

    // sqp_wfqp
    size += sizeof(ocp_nlp_sqp_wfqp_workspace);

    // nlp
    size += ocp_nlp_workspace_calculate_size(config, dims, nlp_opts, in);

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
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    // check for nans
    if (isnan(nlp_res->inf_norm_res_stat) || isnan(nlp_res->inf_norm_res_eq) ||
        isnan(nlp_res->inf_norm_res_ineq) || isnan(nlp_res->inf_norm_res_comp))
    {
        mem->nlp_mem->status = ACADOS_NAN_DETECTED;
        if (nlp_opts->print_level > 0)
        {
            printf("Stopped: NaN detected in iterate.\n");
        }
        return true;
    }

    // check for maximum iterations
    if (!nlp_opts->eval_residual_at_max_iter && n_iter >= nlp_opts->max_iter)
    {
        mem->nlp_mem->status = ACADOS_MAXITER;
        if (nlp_opts->print_level > 0)
        {
            printf("Stopped: Maximum iterations Reached.\n");
        }
        return true;
    }

    // check if solved to tolerance
    if ((nlp_res->inf_norm_res_stat < nlp_opts->tol_stat) &&
        (nlp_res->inf_norm_res_eq < nlp_opts->tol_eq) &&
        (nlp_res->inf_norm_res_ineq < nlp_opts->tol_ineq) &&
        (nlp_res->inf_norm_res_comp < nlp_opts->tol_comp))
    {
        mem->nlp_mem->status = ACADOS_SUCCESS;
        if (nlp_opts->print_level > 0)
        {
            printf("Optimal solution found! Converged to KKT point.\n");
        }
        return true;
    }

    // check for small step
    if (nlp_opts->tol_min_step_norm > 0.0 && (n_iter > 0) && (mem->step_norm < nlp_opts->tol_min_step_norm))
    {
        if (nlp_opts->print_level > 0)
        {
            if (nlp_res->inf_norm_res_eq < nlp_opts->tol_eq && nlp_res->inf_norm_res_ineq < nlp_opts->tol_ineq)
            {
                printf("Stopped: Converged to feasible point. Step size is < tol_eq.\n");
            }
            else
            {
                printf("Stopped: Converged to infeasible point. Step size is < tol_eq.\n");
            }
        }
        mem->nlp_mem->status = ACADOS_MINSTEP;
        return true;
    }

    // check for unbounded problem
    if (mem->nlp_mem->cost_value <= nlp_opts->tol_unbounded)
    {
        mem->nlp_mem->status = ACADOS_UNBOUNDED;
        if (nlp_opts->print_level > 0)
        {
            printf("Stopped: Problem seems to be unbounded.\n");
        }
        return true;
    }

    // check for maximum iterations
    if (n_iter >= nlp_opts->max_iter)
    {
        mem->nlp_mem->status = ACADOS_MAXITER;
        if (nlp_opts->print_level > 0)
        {
            printf("Stopped: Maximum iterations reached.\n");
        }
        return true;
    }

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
    printf("%9s   %9s   %8s   %9s   %9s   ", "step_norm", "step_type", "lm_reg.", "||pi||", "||lam||");
    config->globalization->print_iteration_header();
    printf("\n");
    }
    // print iteration
    ocp_nlp_common_print_iteration(iter, nlp_res);
    printf("%9.2e   %9s   %8.2e   %9.2e   %9.2e   ", mem->step_norm, mem->search_direction_type, prev_levenberg_marquardt, mem->norm_inf_pi, mem->norm_inf_lam);
    config->globalization->print_iteration(nlp_mem->cost_value, nlp_opts->globalization, nlp_mem->globalization);
    printf("\n");
}

/*
Debugging function. Prints:
- ni: number of inequalities, ns: number of slack variables in NLP, nns: number non-slacked constraints (excluding u bounds)
- idxs: indices that were slacked by the user per stage.
- idxns: indices that were slacked by the algorithm for feasibility QP.
*/
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
 * l1 feasibility and pred feasibility functions
 ************************************************/
/*
Calculates predicted reduction of l1 infeasibility by a QP.
- This is defined by l1_inf_QP(0 step) - l1_inf_QP(search direction).
- l1_inf_QP(0 step) == l1_infeasibility at current iterate
- If QP is feasible, then l1_inf_QP(search direction) == 0
*/
static double calculate_pred_l1_inf(ocp_nlp_sqp_wfqp_opts* opts, ocp_nlp_sqp_wfqp_memory *mem, double qp_infeasibility)
{
    if (mem->search_direction_mode == NOMINAL_QP)
    {
        return mem->l1_infeasibility;
    }
    else
    {
        if (mem->l1_infeasibility < MIN(opts->nlp_opts->tol_ineq, opts->nlp_opts->tol_eq))
        {
            return 0.0;
        }
        else
        {
            return mem->l1_infeasibility - qp_infeasibility;
        }
    }
}

/*
Calculates the QP l1 infeasibility by summing up the additional slack variables
included in the QP.
*/
static double calculate_qp_l1_infeasibility_from_slacks(ocp_nlp_dims *dims, ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_sqp_wfqp_opts* opts, ocp_qp_out *qp_out)
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
            l1_inf += MAX(0.0, tmp1);

            // Add upper slack
            tmp2 = BLASFEO_DVECEL(qp_out->ux + i, nx[i]+nu[i]+2*ns[i]+nns[i] + j);
            l1_inf += MAX(0.0, tmp2);
        }
    }
#if defined(ACADOS_DEVELOPER_DEBUG_CHECKS)
    if (l1_inf < -opts->nlp_opts->tol_ineq)
    {
        printf("calculate_qp_l1_infeasibility_from_slacks: got negative l1 infeasibility!\n");
        exit(1);
    }
#endif
    return l1_inf;
}



/*
Calculates the QP l1 infeasibility by explicitely evaluating the QP constraints.
*/
static double calculate_qp_l1_infeasibility_manually(ocp_nlp_dims *dims, ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_sqp_wfqp_workspace *work, ocp_nlp_sqp_wfqp_opts* opts, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
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

        // upper bounds
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
                    l1_inf += MAX(0.0, tmp_bound-tmp);
                    // printf("lower bounds: bound: %.4e, value: %.4e, result: %.4e\n", tmp_bound, tmp, MAX(0.0, tmp_bound-tmp));
                }
                else
                {
                    // upper bounds have the wrong sign!
                    // it is lower_bounds <= value <= -upper_bounds, therefore plus below
                    // printf("upper bounds: value: %.4e, value: %.4e, result: %.4e\n", tmp_bound, tmp, MAX(0.0, tmp_bound+tmp));
                    l1_inf += MAX(0.0, tmp_bound+tmp);
                }
            }
        }
    }
#if defined(ACADOS_DEVELOPER_DEBUG_CHECKS)
    if (l1_inf < -opts->nlp_opts->tol_ineq)
    {
        printf("calculate_qp_l1_infeasibility_manually: got negative l1 infeasibility!\n");
        exit(1);
    }
#endif
    return l1_inf;
}



static double calculate_qp_l1_infeasibility(ocp_nlp_dims *dims, ocp_nlp_sqp_wfqp_memory *mem,
                                            ocp_nlp_sqp_wfqp_workspace *work, ocp_nlp_sqp_wfqp_opts* opts,
                                            ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    double l1_inf = 0.0;
    // this is only possible if directly after a QP was solved
    if (opts->use_QP_l1_inf_from_slacks)
    {
        // Inaccurate if QP solver tolerance is low!
        l1_inf = calculate_qp_l1_infeasibility_from_slacks(dims, mem, opts, qp_out);
    }
    else
    {
        l1_inf = calculate_qp_l1_infeasibility_manually(dims, mem, work, opts, qp_in, qp_out);
    }
    return l1_inf;
}


/************************************************
 * functions for QP preparation
 ************************************************/
/*
Hessian is split into Hessian of cost and Hessian of constraints.
*/
static void set_pointers_for_hessian_evaluation(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
    ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_workspace *work)
{
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    int N = dims->N;

    int *nx = dims->nx;
    int *nu = dims->nu;

    if (!nlp_mem->compute_hess)
    {
        printf("set_pointers_for_hessian_evaluation: constant hessian not supported!\n\n");
        exit(1);
    }

    // init terminal constraint Hessians to 0, as dynamics do not write into them.
    // NOTE: one can implement add_hess_contribution to avoid set
    blasfeo_dgese(nu[N] + nx[N], nu[N] + nx[N], 0.0, mem->RSQ_constr+N, 0, 0);

    // in Hessian computation modules are called in order:
    // dyn, cost, constr.
    int dyn_compute_hess;
    for (int i = 0; i < N; i++)
    {
        config->dynamics[i]->memory_set_RSQrq_ptr(mem->RSQ_constr+i, nlp_mem->dynamics[i]);
        // dynamics always write into hess directly, if hess is computed
        config->dynamics[i]->opts_get(config->dynamics[i], opts->dynamics[i], "compute_hess", &dyn_compute_hess);
        if (!dyn_compute_hess)
        {
            // if dynamics do not compute Hessian, we set it to 0
            blasfeo_dgese(nx[i] + nu[i], nx[i] + nu[i], 0.0, mem->RSQ_constr+i, 0, 0);
        }
    }
    // write cost hess contribution to RSQ_cost
    int add_cost_hess_contribution = 0;
    for (int i = 0; i <= N; i++)
    {
        config->cost[i]->opts_set(config->cost[i], opts->cost[i], "add_hess_contribution", &add_cost_hess_contribution);
        config->cost[i]->memory_set_RSQrq_ptr(mem->RSQ_cost+i, nlp_mem->cost[i]);
    }
    for (int i = 0; i <= N; i++)
    {
        config->constraints[i]->memory_set_RSQrq_ptr(mem->RSQ_constr+i, nlp_mem->constraints[i]);
    }
    return;
}



static void setup_hessian_matrices_for_qps(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_sqp_wfqp_opts *opts,
    ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_workspace *work)
{
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_qp_in *nominal_qp_in = nlp_mem->qp_in;
    ocp_qp_in *relaxed_qp_in = mem->relaxed_qp_in;
    int N = dims->N;

    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;

    int nxu;
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        /* Hessian matrices */
        // hess_QP = hess_cost + hess_constraints
        nxu = nx[i]+nu[i];

        if (opts->use_constraint_hessian_in_feas_qp)
        {
            // Either we use the exact objective Hessian
            blasfeo_dgecp(nxu, nxu, mem->RSQ_constr+i, 0, 0, relaxed_qp_in->RSQrq+i, 0, 0);
        }

        blasfeo_dgecp(nxu, nxu, mem->RSQ_constr+i, 0, 0, nominal_qp_in->RSQrq+i, 0, 0);
        blasfeo_dgead(nxu, nxu, 1.0, mem->RSQ_cost+i, 0, 0, nominal_qp_in->RSQrq+i, 0, 0);

        // Z -- slack matrix --> needs to be at correct position!
        blasfeo_dveccpsc(2*ns[i], 1.0, mem->Z_cost_module+i, 0, nominal_qp_in->Z+i, 0);
    }
    // Levenberg Marquardt term for nominal QP
    ocp_nlp_add_levenberg_marquardt_term(config, dims, in, out, opts->nlp_opts, nlp_mem, work, mem->alpha, nlp_mem->iter, nlp_mem->qp_in);
}

/*
Solves the QP. Either solves feasibility QP or nominal QP
*/
static int prepare_and_solve_QP(ocp_nlp_config* config, ocp_nlp_sqp_wfqp_opts* opts,
                    ocp_qp_in* scaled_qp_in, ocp_qp_in* qp_in, ocp_qp_out* scaled_qp_out, ocp_qp_out* qp_out,
                    ocp_nlp_dims *dims, ocp_nlp_sqp_wfqp_memory* mem, ocp_nlp_in* nlp_in, ocp_nlp_out* nlp_out,
                    ocp_nlp_memory* nlp_mem, ocp_nlp_workspace* nlp_work, bool solve_feasibility_qp,
                    acados_timer timer_tot)
{
    acados_timer timer;
    ocp_nlp_opts* nlp_opts = opts->nlp_opts;
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_timings *nlp_timings = nlp_mem->nlp_timings;
    int qp_status = ACADOS_SUCCESS;

    // printf("\nprepare_and_solve_QP: solve_feasibility_qp %d\n", solve_feasibility_qp);
    // warm start of first QP
    if (nlp_mem->iter == 0)
    {
        if (!nlp_opts->warm_start_first_qp)
        {
            // (typically) no warm start at first iteration
            int tmp_int = 0;
            qp_solver->opts_set(qp_solver, nlp_opts->qp_solver_opts, "warm_start", &tmp_int);
        }
        else if (nlp_opts->warm_start_first_qp_from_nlp)
        {
            int tmp_bool = true;
            qp_solver->opts_set(qp_solver, nlp_opts->qp_solver_opts, "initialize_next_xcond_qp_from_qp_out", &tmp_bool);
            ocp_nlp_initialize_qp_from_nlp(config, dims, qp_in, nlp_out, qp_out);
        }
    }

    if (mem->qps_solved_in_iter < 2 && (!solve_feasibility_qp || opts->use_constraint_hessian_in_feas_qp))
    {
        if (solve_feasibility_qp && opts->use_constraint_hessian_in_feas_qp)
        {
            // LM for feasibility QP
            ocp_nlp_add_levenberg_marquardt_term(config, dims, nlp_in, nlp_out, opts->nlp_opts, nlp_mem, nlp_work, mem->alpha, nlp_mem->iter, qp_in);
        }

        acados_tic(&timer);
        // regularize Hessian
        config->regularize->regularize(config->regularize, dims->regularize, nlp_opts->regularize, nlp_mem->regularize_mem);
        nlp_timings->time_reg += acados_toc(&timer);
    }

    // Show input to QP
    if (nlp_opts->print_level > 3)
    {
        printf("\n\nSQP: ocp_qp_in at iteration %d\n", nlp_mem->iter);
        print_ocp_qp_dims(qp_in->dim);
        print_ocp_qp_in(scaled_qp_in);
    }

#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
    ocp_nlp_dump_qp_in_to_file(qp_in, nlp_mem->iter, 0);
#endif
    if (solve_feasibility_qp)
    {
        if (opts->use_constraint_hessian_in_feas_qp)
        {
            qp_status = ocp_nlp_solve_qp_and_correct_dual(config, dims, nlp_opts,
                nlp_mem, nlp_work, false,
                scaled_qp_in, qp_in, scaled_qp_out, qp_out, &mem->relaxed_qp_solver);
        }
        else
        {
            // dont regularize Hessian for feasibility QP
            qp_status = ocp_nlp_solve_qp(config, dims, nlp_opts,
                nlp_mem, nlp_work, scaled_qp_in, scaled_qp_out, &mem->relaxed_qp_solver);
            acados_tic(&timer);
            ocp_nlp_qpscaling_rescale_solution(dims->relaxed_qpscaling, nlp_opts->qpscaling, mem->relaxed_qpscaling_mem, qp_in, qp_out);
            nlp_timings->time_qpscaling += acados_toc(&timer);
        }
    }
    else
    {
        qp_status = ocp_nlp_solve_qp_and_correct_dual(config, dims, nlp_opts,
                                                    nlp_mem, nlp_work, false,
                                                    NULL, NULL, NULL, NULL, NULL);
    }
    mem->qps_solved_in_iter += 1;

    // restore default warm start
    if (nlp_mem->iter==0)
    {
        qp_solver->opts_set(qp_solver, nlp_opts->qp_solver_opts, "warm_start", &nlp_opts->qp_warm_start);
    }

    if (nlp_opts->print_level > 3)
    {
        printf("\n\nSQP: ocp_qp_out at iteration %d\n", nlp_mem->iter);
        print_ocp_qp_dims(qp_out->dim);
        print_ocp_qp_out(scaled_qp_out);
    }

#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
    ocp_nlp_dump_qp_out_to_file(qp_out, nlp_mem->iter, 0);
#endif

    // exit conditions on QP status
    if (qp_status!=ACADOS_SUCCESS)
    {
        if (nlp_opts->print_level > 1)
        {
            printf("\n Failed to solve the following QP:\n");
            if (nlp_opts->print_level)
                print_ocp_qp_in(qp_in);
        }

        mem->nlp_mem->status = ACADOS_QP_FAILURE;
        nlp_timings->time_tot = acados_toc(&timer_tot);

        return mem->nlp_mem->status;
    }
    return qp_status;
}

static void log_qp_stats(ocp_nlp_sqp_wfqp_memory *mem, bool solve_feasibility_qp,
    int qp_status, int qp_iter)
{
    int nlp_iter = mem->nlp_mem->iter;
    if (nlp_iter < mem->stat_m)
    {
        if (mem->search_direction_mode == NOMINAL_QP)
        {
            mem->stat[mem->stat_n*(nlp_iter)+4] = qp_status;
            mem->stat[mem->stat_n*(nlp_iter)+5] = qp_iter;
        }
        else if (mem->search_direction_mode == BYRD_OMOJOKUN)
        {
            if (solve_feasibility_qp)
            {
                mem->stat[mem->stat_n*(nlp_iter)+6] = qp_status;
                mem->stat[mem->stat_n*(nlp_iter)+7] = qp_iter;
            }
            else
            {
                mem->stat[mem->stat_n*(nlp_iter)+8] = qp_status;
                mem->stat[mem->stat_n*(nlp_iter)+9] = qp_iter;
            }
        }
    }
}

static void log_multiplier_norms(ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_sqp_wfqp_opts *opts, ocp_nlp_out *nlp_out, ocp_nlp_dims *dims)
{
    int nlp_iter = mem->nlp_mem->iter;
    if (opts->log_pi_norm_inf)
    {
        mem->norm_inf_pi = ocp_nlp_compute_dual_pi_norm_inf(dims, nlp_out);
        mem->stat[mem->stat_n*(nlp_iter)+11] = mem->norm_inf_pi;
    }
    if (opts->log_lam_norm_inf)
    {
        mem->norm_inf_lam = ocp_nlp_compute_dual_lam_norm_inf(dims, nlp_out);
        mem->stat[mem->stat_n*(nlp_iter)+12] = mem->norm_inf_lam;
    }
}

/************************************************
* Byrd-Omojokun Subproblem Functions:
************************************************/
/*
Adjusts the bounds of the nominal QP with the optimal slack variables of the
feasibility. Resulting nominal QP has always a feasible solution.
*/
static void setup_byrd_omojokun_bounds(ocp_nlp_dims *dims, ocp_nlp_memory *nlp_mem, ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_sqp_wfqp_workspace *work, ocp_nlp_sqp_wfqp_opts* opts)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int *nns = mem->nns;

    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ni_nl = dims->ni_nl;

    ocp_qp_in *nominal_qp_in = nlp_mem->scaled_qp_in;
    ocp_qp_out *relaxed_qp_out = mem->relaxed_scaled_qp_out;

    int i, j;
    double tmp_lower, tmp_upper;

    for (i = 0; i <= N; i++)
    {
        for (j=0; j<mem->nns[i]; ++j)
        {
            int constr_index = mem->idxns[i][j]; //
            int slack_index = mem->relaxed_qp_in->idxs_rev[i][constr_index];
            // get lower slack
            tmp_lower = BLASFEO_DVECEL(relaxed_qp_out->ux + i, nx[i]+nu[i]+slack_index);
            // lower_bound - value
            BLASFEO_DVECEL(nominal_qp_in->d+i, constr_index) -= tmp_lower;

            // get upper slack
            tmp_upper = BLASFEO_DVECEL(relaxed_qp_out->ux + i, nx[i]+nu[i]+ns[i]+nns[i] + slack_index);
            // upper bounds have the wrong sign!
            // it is lower_bounds <= value <= -upper_bounds, therefore plus below
            // for the slacks with upper bound we have value - slack, therefore
            // value <= -upper_bound + slack,
            // we store upper_bound - slack
            BLASFEO_DVECEL(nominal_qp_in->d+i, nb[i] + ng[i] + ni_nl[i] + constr_index) -= tmp_upper;
        }
    }
}

/*
First: solve feasibility QP. Second: solve nominal QP with adjusted bounds.

Guarantees well-defined search direction in case of infeasible nominal QP.
*/
static int byrd_omojokun_direction_computation(ocp_nlp_dims *dims,
                                            ocp_nlp_config *config,
                                            ocp_nlp_sqp_wfqp_opts *opts,
                                            ocp_nlp_opts *nlp_opts,
                                            ocp_nlp_in *nlp_in,
                                            ocp_nlp_out *nlp_out,
                                            ocp_nlp_sqp_wfqp_memory *mem,
                                            ocp_nlp_sqp_wfqp_workspace *work,
                                            acados_timer timer_tot)
{
    ocp_nlp_memory* nlp_mem = mem->nlp_mem;
    ocp_nlp_workspace* nlp_work = work->nlp_work;

    ocp_qp_in *nominal_qp_in = nlp_mem->qp_in;
    ocp_qp_in *nominal_scaled_qp_in = nlp_mem->scaled_qp_in;
    ocp_qp_out *nominal_qp_out = nlp_mem->qp_out;
    ocp_qp_out *nominal_scaled_qp_out = nlp_mem->scaled_qp_out;

    ocp_qp_in *relaxed_qp_in = mem->relaxed_qp_in;
    ocp_qp_in *relaxed_scaled_qp_in = mem->relaxed_scaled_qp_in;
    ocp_qp_out *relaxed_qp_out = mem->relaxed_qp_out;

    ocp_nlp_timings *nlp_timings = nlp_mem->nlp_timings;
    qp_info* qp_info_;

    int qp_status;
    int qp_iter = 0;

    /* Solve Feasibility QP: Objective: Only constraint Hessian/Identity AND only gradient of slack variables */
    print_debug_output("Solve Feasibility QP!\n", nlp_opts->print_level, 2);
    qp_status = prepare_and_solve_QP(config, opts, relaxed_scaled_qp_in, relaxed_qp_in, mem->relaxed_scaled_qp_out, relaxed_qp_out, dims, mem, nlp_in, nlp_out,
                nlp_mem, nlp_work, true, timer_tot);
    ocp_qp_out_get(relaxed_qp_out, "qp_info", &qp_info_);
    qp_iter = qp_info_->num_iter;
    log_qp_stats(mem, true, qp_status, qp_iter);
    if (qp_status != ACADOS_SUCCESS)
    {
        if (nlp_opts->print_level >=1)
        {
            printf("\nError in feasibility QP in iteration %d, got qp_status %d!\n", qp_iter, qp_status);
        }
        nlp_mem->status = ACADOS_QP_FAILURE;
        nlp_timings->time_tot = acados_toc(&timer_tot);
        return nlp_mem->status;
    }

    // here was the calculation of pred_infeasibility earlier, moved outside

    /* Solve the nominal QP with updated bounds*/
    print_debug_output("Solve Nominal QP!\n", nlp_opts->print_level, 2);
    setup_byrd_omojokun_bounds(dims, nlp_mem, mem, work, opts);
    // solve_feasibility_qp --> false in prepare_and_solve_QP

    qp_status = prepare_and_solve_QP(config, opts, nominal_scaled_qp_in, nominal_qp_in, nominal_scaled_qp_out, nominal_qp_out, dims, mem, nlp_in, nlp_out,
                                     nlp_mem, nlp_work, false, timer_tot);
    ocp_qp_out_get(nominal_qp_out, "qp_info", &qp_info_);
    qp_iter = qp_info_->num_iter;
    log_qp_stats(mem, false, qp_status, qp_iter);

    if (qp_status != ACADOS_SUCCESS)
    {
        if (nlp_opts->print_level >=1)
        {
            printf("\nError in nominal QP in iteration %d, got qp_status %d!\n", qp_iter, qp_status);
        }
        nlp_mem->status = ACADOS_QP_FAILURE;
        nlp_timings->time_tot = acados_toc(&timer_tot);
        return nlp_mem->status;
    }
    return qp_status;
}

/********************************
* Functions for relaxed QP:

The relaxed QP is designed such that it is always feasible. We achieve this by
adding slack variables to all constraints that were not slacked by the user upfront.
Exceptations to that slacking are:
- bounds on states at the initial node
- all bounds on controls
- dynamics (since always full rank)

Hessian matrix is either identity for states or the hessian of constraints given by
the Hessian mode
Feasibility QP only has the gradient of the slack variables, but NOT the gradient of the objective.

Aim of feasibility QP: Achieve maximum reduction of infeasibility of QP, if nominal QP is not feasible

dynamics: shared with nominal QP, set_relaxed_qp_in_matrix_pointers
- BAbt
- b

constraint definitions: shared with nominal QP, set_relaxed_qp_in_matrix_pointers
- DCt
- idxb
- idxe


- d_mask: different due to slacks in feasibility QP setup once in ocp_nlp_sqp_wfqp_approximate_feasibility_qp_constraint_vectors

Hessian:
- RSQrq: if identity in initial_setup_feasibility_qp_objective, otherwise in setup_hessian_matrices_for_qps
- Z: set in initial_setup_feasibility_qp_objective

Vectors
- rqz: initial_setup_feasibility_qp_objective
- d: ocp_nlp_sqp_wfqp_approximate_feasibility_qp_constraint_vectors

Index vectors
- idxs_rev: set in ocp_nlp_sqp_wfqp_precompute
- diag_H_flag : set in initial_setup_feasibility_qp_objective

- m not set at the moment

*********************************/
/*
Feasibility QP and nominal QP share many entries.
- Constraint matrices are always the same
- Cost gradient and Hessian are different
- Constraint bounds are different
Where we point to the same memory for both QPs is given below.
*/
void set_relaxed_qp_in_matrix_pointers(ocp_nlp_sqp_wfqp_memory *mem)
{
    ocp_qp_in *qp_in = mem->nlp_mem->qp_in;

    // dynamics
    mem->relaxed_qp_in->BAbt = qp_in->BAbt; // dynamics matrix & vector work space
    mem->relaxed_qp_in->b = qp_in->b; // dynamics vector

    // constraint defintitions
    mem->relaxed_qp_in->DCt = qp_in->DCt; // inequality constraints matrix
    mem->relaxed_qp_in->idxb = qp_in->idxb;
    mem->relaxed_qp_in->idxe = qp_in->idxe;
}

static void set_relaxed_scaled_qp_in_matrix_pointers(ocp_nlp_sqp_wfqp_memory *mem)
{
    ocp_qp_in *scaled_qp_in = mem->nlp_mem->scaled_qp_in;

    // dynamics
    mem->relaxed_scaled_qp_in->BAbt = scaled_qp_in->BAbt;
    mem->relaxed_scaled_qp_in->b = scaled_qp_in->b;
    // constraint defintitions
    mem->relaxed_scaled_qp_in->DCt = scaled_qp_in->DCt;
    mem->relaxed_scaled_qp_in->idxb = scaled_qp_in->idxb;
    mem->relaxed_scaled_qp_in->idxe = scaled_qp_in->idxe;
}

/*
update QP rhs for feasibility QP (step prim var, abs dual var)
- copy d parts from nominal QP and include bounds of QP slacks
- setup d_mask for QP slacks
*/
void ocp_nlp_sqp_wfqp_approximate_feasibility_qp_constraint_vectors(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
    ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_workspace *work, bool scaled)
{
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_qp_in *nominal_qp_in;
    ocp_qp_in *relaxed_qp_in;
    if (scaled)
    {
        nominal_qp_in = nlp_mem->scaled_qp_in;
        relaxed_qp_in = mem->relaxed_scaled_qp_in;
    }
    else
    {
        nominal_qp_in = nlp_mem->qp_in;
        relaxed_qp_in = mem->relaxed_qp_in;
    }
    int N = dims->N;

    int *ns = dims->ns;
    int *ni = dims->ni;
    int *nns = mem->nns;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // d --> copy what is possible from nominal_qp_in
        int n_nominal_ineq_nlp = ni[i] - ns[i];
        blasfeo_dveccp(2*n_nominal_ineq_nlp+ns[i], nominal_qp_in->d + i, 0, relaxed_qp_in->d + i, 0);
        blasfeo_dveccp(ns[i], nominal_qp_in->d + i, 2*n_nominal_ineq_nlp+ns[i], relaxed_qp_in->d + i, 2*n_nominal_ineq_nlp+ns[i]+nns[i]);
    }
    // setup d_mask
    if (nlp_mem->iter == 0)
    {
        int offset_dmask;
        for (int i=0; i<=dims->N; i++)
        {
            offset_dmask = 2*(dims->nb[i]+dims->ng[i]+dims->ni_nl[i]);
            blasfeo_dveccp(offset_dmask, nominal_qp_in->d_mask+i, 0, relaxed_qp_in->d_mask+i, 0);
            blasfeo_dvecse(2*relaxed_qp_in->dim->ns[i], 1.0, relaxed_qp_in->d_mask+i, offset_dmask);
        }
    }
}


/*
Sets Hessian and gradient of feasibility QP at start of solve process.
- gradient is always constant for all feasibility QPs
- If identity Hessian is used, then Hessian in feasibility QP is also constant
*/
static void initial_setup_feasibility_qp_objective(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_sqp_wfqp_opts *opts,
    ocp_nlp_sqp_wfqp_memory *mem, ocp_nlp_workspace *work)
{
    ocp_qp_in *relaxed_qp_in = mem->relaxed_qp_in;

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
        // hess_QP = hess_cost + hess_constraints
        nxu = nx[i]+nu[i];

        blasfeo_dgese(nxu, nxu, 0.0, relaxed_qp_in->RSQrq+i, 0, 0);
        if (!opts->use_constraint_hessian_in_feas_qp)
        {
            // We use the identity matrix Hessian
            blasfeo_ddiare(nxu, opts->feasibility_qp_hessian_scalar, relaxed_qp_in->RSQrq+i, 0, 0);
            // indicate that Hessian is diagonal
            relaxed_qp_in->diag_H_flag[i] = 1;
        }

        // Z -- slack matrix; order: Z = [Zl_NLP, Zl_QP, Zu_NLP, Zu_QP]

        // Zl_QP
        blasfeo_dvecse(nns[i], 0.0, relaxed_qp_in->Z+i, ns[i]);
        // Zu_QP
        blasfeo_dvecse(nns[i], 0.0, relaxed_qp_in->Z+i, 2*ns[i]+nns[i]);

        if (opts->use_constraint_hessian_in_feas_qp)
        {
            // Zl_NLP
            blasfeo_dvecse(ns[i], 0.0, relaxed_qp_in->Z+i, 0);
            // Zu_NLP
            blasfeo_dvecse(ns[i], 0.0, relaxed_qp_in->Z+i, ns[i]+nns[i]);
        }
        else
        {
            // Zl_NLP
            blasfeo_dvecse(ns[i], opts->feasibility_qp_hessian_scalar, relaxed_qp_in->Z+i, 0);
            // Zu_NLP
            blasfeo_dvecse(ns[i], opts->feasibility_qp_hessian_scalar, relaxed_qp_in->Z+i, ns[i]+nns[i]);
        }

        /* vectors */
        // rqz
        // be aware of rqz_QP = [r, q, zl_NLP, zl_QP, zu_NLP, zu_QP]
        blasfeo_dvecse(nx[i]+nu[i]+ns[i], 0.0, relaxed_qp_in->rqz+i, 0);
        blasfeo_dvecse(ns[i], 0.0, relaxed_qp_in->rqz+i, nx[i]+nu[i]+ns[i]+nns[i]);

        // zl_QP
        blasfeo_dvecse(nns[i], 1.0, relaxed_qp_in->rqz+i, nu[i]+nx[i]+ns[i]);
        // zu_QP
        blasfeo_dvecse(nns[i], 1.0, relaxed_qp_in->rqz+i, nu[i]+nx[i]+2*ns[i]+nns[i]);

    }
}

/************************************************
* Search Direction Function
************************************************/
/*
Depending on the search direction mode, a search direction is calculated.

- NOMINAL_QP solves nominal QP (might be infeasible)
- if NOMINAL_QP infeasible, switch to BYRD_OMOJOKUN mode
- BYRD_OMOJOKUN guarantees always well-defined search directions
- If slack variables are 0 in feasibility for enough consecutive iterations switch back to NOMINAL_QP mode
*/
static int calculate_search_direction(ocp_nlp_dims *dims,
    ocp_nlp_config *config,
    ocp_nlp_sqp_wfqp_opts *opts,
    ocp_nlp_opts *nlp_opts,
    ocp_nlp_in *nlp_in,
    ocp_nlp_out *nlp_out,
    ocp_nlp_sqp_wfqp_memory *mem,
    ocp_nlp_sqp_wfqp_workspace *work,
    acados_timer timer_tot)
{
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    qp_info* qp_info_;
    int qp_iter = 0;
    int search_direction_status;
    mem->qps_solved_in_iter = 0;

    if (mem->search_direction_mode == NOMINAL_QP)
    {
        // if the QP can be solved and the status is good, we return 0
        // otherwise, we change the mode to Byrd-Omojokun and we continue.
        search_direction_status = prepare_and_solve_QP(config, opts, nlp_mem->scaled_qp_in,
            nlp_mem->qp_in, nlp_mem->scaled_qp_out, nlp_mem->qp_out,
            dims, mem, nlp_in, nlp_out, nlp_mem, work->nlp_work, false, timer_tot);
        ocp_qp_out_get(nlp_mem->qp_out, "qp_info", &qp_info_);
        qp_iter = qp_info_->num_iter;
        log_qp_stats(mem, false, search_direction_status, qp_iter);
        if (search_direction_status != ACADOS_SUCCESS)
        {
            if (nlp_opts->print_level >1)
            {
                printf("\nError in nominal QP in iteration %d, got qp_status %d!\n", qp_iter, search_direction_status);
                printf("Switch to Byrd-Omojokun mode!\n");
            }
            mem->search_direction_mode = BYRD_OMOJOKUN;
        }
        else
        {
            mem->search_direction_type = "N";
            if (config->globalization->needs_objective_value() == 1)
            {
                /*
                Calculates predicted reduction of l1 infeasibility by a QP. This is defined by
                l1_inf_QP(0 step) - search direction)
                */
                mem->pred_l1_inf_QP = calculate_pred_l1_inf(opts, mem, -1.0);
            }
            return ACADOS_SUCCESS;
        }
    }
    if (mem->search_direction_mode == BYRD_OMOJOKUN)
    {
        // We solve two QPs and return the search direction that we found!
        // if the second QP is feasible, we change back to nominal QP mode.
        if (mem->qps_solved_in_iter == 1)
        {
            mem->search_direction_type = "NFN";
        }
        else
        {
            mem->search_direction_type = "FN";
        }
        search_direction_status = byrd_omojokun_direction_computation(dims, config, opts, nlp_opts, nlp_in, nlp_out, mem, work, timer_tot);

        double l1_inf_QP_feasibility = calculate_qp_l1_infeasibility(dims, mem, work, opts, mem->relaxed_qp_in, mem->relaxed_qp_out);
        if (config->globalization->needs_objective_value() == 1)
        {
            mem->pred_l1_inf_QP = calculate_pred_l1_inf(opts, mem, l1_inf_QP_feasibility);
        }

        if (l1_inf_QP_feasibility/(MAX(1.0, (double) mem->absolute_nns)) < nlp_opts->tol_ineq)
        {
            mem->watchdog_zero_slacks_counter += 1;
        }

        if (opts->allow_direction_mode_switch_to_nominal && mem->watchdog_zero_slacks_counter == opts->watchdog_zero_slacks_max)
        {
            mem->watchdog_zero_slacks_counter = 0;
            mem->search_direction_mode = NOMINAL_QP;
        }

        return search_direction_status;
        // Maybe switch to a full feasibility restoration phase if NLP seems infeasible? Will be implemented later
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

/************************************************
* MAIN OPTIMIZATION ROUTINE
************************************************/

int ocp_nlp_sqp_wfqp(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_)
{
    acados_timer timer_tot, timer1;
    acados_tic(&timer_tot);

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

    // zero timers
    ocp_nlp_timings_reset(nlp_timings);

    int qp_status = 0;
    int qp_iter = 0;
    nlp_mem->objective_multiplier = 1.0;
    mem->alpha = 0.0;
    mem->step_norm = 0.0;
    mem->norm_inf_lam = 0.0;
    mem->norm_inf_pi = 0.0;
    mem->nlp_mem->status = ACADOS_READY;
    mem->search_direction_type = "-";
    mem->search_direction_mode = opts->search_direction_mode;
    mem->watchdog_zero_slacks_counter = 0;
    mem->l1_infeasibility = -1.0; // default, cannot be negative

#if defined(ACADOS_WITH_OPENMP)
    // backup number of threads
    int num_threads_bkp = omp_get_num_threads();
    // set number of threads
    omp_set_num_threads(opts->nlp_opts->num_threads);
#endif

    ocp_nlp_initialize_submodules(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);

    // gradient of feasibility QP is always constant. So is Hessian, if identity Hessian is used
    initial_setup_feasibility_qp_objective(config, dims, nlp_in, nlp_out, opts, mem, nlp_work);

    /************************************************
     * main sqp loop
     ************************************************/
    nlp_mem->iter = 0;
    double prev_levenberg_marquardt = 0.0;
    int search_direction_status = 0;

    if (nlp_opts->print_level > 1)
    {
        print_indices(dims, work, mem);
    }

    for (; nlp_mem->iter <= opts->nlp_opts->max_iter; nlp_mem->iter++) // <= needed such that after last iteration KKT residuals are checked before max_iter is thrown.
    {
        // We always evaluate the residuals until the last iteration
        // If the option "eval_residual_at_max_iter" is set, we also
        // evaluate the residuals after the last iteration.
        if (nlp_mem->iter != opts->nlp_opts->max_iter || nlp_opts->eval_residual_at_max_iter)
        {
            // store current iterate
            if (nlp_opts->store_iterates)
            {
                copy_ocp_nlp_out(dims, nlp_out, nlp_mem->iterates[nlp_mem->iter]);
            }
            /* Prepare the QP data */
            // linearize NLP and update QP matrices
            acados_tic(&timer1);
            set_pointers_for_hessian_evaluation(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work);
            // nominal QP solver
            ocp_nlp_approximate_qp_matrices(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
            ocp_nlp_approximate_qp_vectors_sqp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);

            if (nlp_opts->with_adaptive_levenberg_marquardt || config->globalization->needs_objective_value() == 1)
            {
                ocp_nlp_get_cost_value_from_submodules(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
            }

            // setup relaxed QP
            // matrices for relaxed QP solver evaluated in nominal QP solver
            ocp_nlp_sqp_wfqp_approximate_feasibility_qp_constraint_vectors(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work, false);
            setup_hessian_matrices_for_qps(config, dims, nlp_in, nlp_out, opts, mem, nlp_work);
            //
            nlp_timings->time_lin += acados_toc(&timer1);

            // compute nlp residuals
            ocp_nlp_res_compute(dims, nlp_opts, nlp_in, nlp_out, nlp_res, nlp_mem, nlp_work);
            ocp_nlp_res_get_inf_norm(nlp_res, &nlp_out->inf_norm_res);
        }

        // Initialize the memory for different globalization strategies
        if (nlp_mem->iter == 0)
        {
            config->globalization->initialize_memory(config, dims, nlp_mem, nlp_opts);
        }

        // save statistics
        if (nlp_mem->iter < mem->stat_m)
        {
            mem->stat[mem->stat_n*nlp_mem->iter+0] = nlp_res->inf_norm_res_stat;
            mem->stat[mem->stat_n*nlp_mem->iter+1] = nlp_res->inf_norm_res_eq;
            mem->stat[mem->stat_n*nlp_mem->iter+2] = nlp_res->inf_norm_res_ineq;
            mem->stat[mem->stat_n*nlp_mem->iter+3] = nlp_res->inf_norm_res_comp;
        }

        log_multiplier_norms(mem, opts, nlp_out, dims);

        /* Output */
        if (nlp_opts->print_level > 0)
        {
            print_iteration(nlp_mem->iter, config, nlp_res, mem, nlp_opts, prev_levenberg_marquardt, qp_status, qp_iter);
        }
        prev_levenberg_marquardt = nlp_opts->levenberg_marquardt;

        /* Termination */
        if (check_termination(nlp_mem->iter, dims, nlp_res, mem, opts))
        {
#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif
            nlp_timings->time_tot = acados_toc(&timer_tot);
            return mem->nlp_mem->status;
        }

        if (config->globalization->needs_objective_value() == 1)
        {
            mem->l1_infeasibility = ocp_nlp_get_l1_infeasibility(config, dims, nlp_mem);
        }

        /* Scale the QP */
        // scale the qp: includes constraints and objective
        acados_tic(&timer1);
        ocp_nlp_qpscaling_scale_qp(dims->qpscaling, nlp_opts->qpscaling, nlp_mem->qpscaling, nominal_qp_in);
        nlp_timings->time_qpscaling += acados_toc(&timer1);
        ocp_nlp_sqp_wfqp_approximate_feasibility_qp_constraint_vectors(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work, true);

        acados_tic(&timer1);
        ocp_nlp_qpscaling_scale_qp(dims->relaxed_qpscaling, nlp_opts->qpscaling, mem->relaxed_qpscaling_mem, mem->relaxed_qp_in); // ensures feasibility constraint Hessian is scaled
        nlp_timings->time_qpscaling += acados_toc(&timer1);

        /* Search Direction Computation */
        search_direction_status = calculate_search_direction(dims, config, opts, nlp_opts, nlp_in, nlp_out, mem, work, timer_tot);
        if (search_direction_status != ACADOS_SUCCESS)
        {
#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif
            return nlp_mem->status;
        }

        // Compute the step norm
        if (nlp_opts->tol_min_step_norm > 0.0 || nlp_opts->log_primal_step_norm)
        {
            mem->step_norm = ocp_qp_out_compute_primal_nrm_inf(nominal_qp_out);
            if (nlp_opts->log_primal_step_norm)
                nlp_mem->primal_step_norm[nlp_mem->iter] = mem->step_norm;
        }
        if (nlp_opts->log_dual_step_norm && !nlp_opts->with_anderson_acceleration)
        {
            nlp_mem->dual_step_norm[nlp_mem->iter] = ocp_nlp_compute_delta_dual_norm_inf(dims, nlp_work, nlp_out, nominal_qp_out);
        }

        /* globalization */
        // Calculate optimal QP objective (needed for globalization)
        if (config->globalization->needs_qp_objective_value() == 1)
        {
            nlp_mem->qp_cost_value = ocp_nlp_compute_qp_objective_value(dims, nominal_qp_in, nominal_qp_out, nlp_work);
            nlp_mem->predicted_infeasibility_reduction = mem->pred_l1_inf_QP;
            nlp_mem->predicted_optimality_reduction = -ocp_nlp_compute_gradient_directional_derivative(dims, nominal_qp_in, nominal_qp_out);
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
            nlp_timings->time_tot = acados_toc(&timer_tot);
#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif
            return nlp_mem->status;
        }

        mem->stat[mem->stat_n*(nlp_mem->iter+1)+10] = mem->alpha;
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
    ocp_qp_xcond_solver relaxed_qp_solver = mem->relaxed_qp_solver;

    config->qp_solver->memory_reset(qp_solver, dims->qp_solver,
        nlp_mem->qp_in, nlp_mem->qp_out, opts->nlp_opts->qp_solver_opts,
        nlp_mem->qp_solver_mem, nlp_work->qp_work);
    relaxed_qp_solver.config->memory_reset(relaxed_qp_solver.config,
                                           relaxed_qp_solver.dims,
                                           mem->relaxed_qp_in,
                                           mem->relaxed_qp_out,
                                           opts->nlp_opts->qp_solver_opts,
                                           mem->relaxed_qp_solver_mem,
                                           mem->relaxed_qp_solver_work);
}


int ocp_nlp_sqp_wfqp_setup_qp_matrices_and_factorize(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
    void *opts_, void *mem_, void *work_)
{
    ocp_nlp_sqp_wfqp_opts *opts = opts_;
    ocp_nlp_sqp_wfqp_memory *mem = mem_;
    ocp_nlp_sqp_wfqp_workspace *work = work_;

    return ocp_nlp_common_setup_qp_matrices_and_factorize(config_, dims_, nlp_in_, nlp_out_, opts->nlp_opts, mem->nlp_mem, work->nlp_work);
}


void ocp_nlp_sqp_wfqp_eval_kkt_residual(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
    void *opts_, void *mem_, void *work_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_sqp_wfqp_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_nlp_sqp_wfqp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    ocp_nlp_initialize_submodules(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
    ocp_nlp_approximate_qp_matrices(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
    ocp_nlp_approximate_qp_vectors_sqp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
    ocp_nlp_res_compute(dims, nlp_opts, nlp_in, nlp_out, nlp_mem->nlp_res, nlp_mem, nlp_work);
}




/*
Computes indices of constraints in NLP that were not slacked by the user. Excludes
bounds on controls (u)
*/
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

    ocp_nlp_qpscaling_precompute(dims->relaxed_qpscaling, opts->nlp_opts->qpscaling, mem->relaxed_qpscaling_mem, mem->relaxed_qp_in, mem->relaxed_qp_out);
    ocp_nlp_qpscaling_memory_get(dims->relaxed_qpscaling, mem->relaxed_qpscaling_mem, "scaled_qp_in", 0, &mem->relaxed_scaled_qp_in);
    ocp_nlp_qpscaling_memory_get(dims->relaxed_qpscaling, mem->relaxed_qpscaling_mem, "scaled_qp_out", 0, &mem->relaxed_scaled_qp_out);

    set_relaxed_scaled_qp_in_matrix_pointers(mem);

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
    acados_timer timer_tot;
    acados_tic(&timer_tot);

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

    nlp_mem->nlp_timings->time_solution_sensitivities = acados_toc(&timer_tot);

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

void ocp_nlp_sqp_wfqp_eval_solution_sens_adj_p(void *config_, void *dims_,
    void *opts_, void *mem_, void *work_, void *sens_nlp_out,
    const char *field, int stage, void *grad_p)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_wfqp_opts *opts = opts_;
    ocp_nlp_sqp_wfqp_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_nlp_sqp_wfqp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;
    ocp_nlp_common_eval_solution_sens_adj_p(config, dims,
        opts->nlp_opts, nlp_mem, nlp_work,
        sens_nlp_out, field, stage, grad_p);
}


void ocp_nlp_sqp_wfqp_get(void *config_, void *dims_, void *mem_, const char *field, void *return_value_)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_sqp_wfqp_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

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
        ocp_nlp_timings_get(config, nlp_mem->nlp_timings, field, return_value_);
    }
    else if (!strcmp("stat", field))
    {
        double **value = return_value_;
        *value = mem->stat;
    }
    else if (!strcmp("statistics", field))
    {
        int n_row = mem->stat_m<nlp_mem->iter+1 ? mem->stat_m : nlp_mem->iter+1;
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
        ocp_nlp_memory_get(config, nlp_mem, field, return_value_);
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
    config->setup_qp_matrices_and_factorize = &ocp_nlp_sqp_wfqp_setup_qp_matrices_and_factorize;
    config->memory_reset_qp_solver = &ocp_nlp_sqp_wfqp_memory_reset_qp_solver;
    config->eval_param_sens = &ocp_nlp_sqp_wfqp_eval_param_sens;
    config->eval_lagr_grad_p = &ocp_nlp_sqp_wfqp_eval_lagr_grad_p;
    config->eval_solution_sens_adj_p = &ocp_nlp_sqp_wfqp_eval_solution_sens_adj_p;
    config->config_initialize_default = &ocp_nlp_sqp_wfqp_config_initialize_default;
    config->precompute = &ocp_nlp_sqp_wfqp_precompute;
    config->get = &ocp_nlp_sqp_wfqp_get;
    config->opts_get = &ocp_nlp_sqp_wfqp_opts_get;
    config->work_get = &ocp_nlp_sqp_wfqp_work_get;
    config->terminate = &ocp_nlp_sqp_wfqp_terminate;
    config->step_update = &ocp_nlp_update_variables_sqp;
    config->is_real_time_algorithm = &ocp_nlp_sqp_wfqp_is_real_time_algorithm;
    config->eval_kkt_residual = &ocp_nlp_sqp_wfqp_eval_kkt_residual;

    return;
}
