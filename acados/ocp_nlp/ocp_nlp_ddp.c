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


#include "acados/ocp_nlp/ocp_nlp_ddp.h"

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
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados_c/ocp_qp_interface.h"



/************************************************
 * options
 ************************************************/

acados_size_t ocp_nlp_ddp_opts_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_ddp_opts);

    size += ocp_nlp_opts_calculate_size(config, dims);

    return size;
}



void *ocp_nlp_ddp_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_ddp_opts *opts = (ocp_nlp_ddp_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_ddp_opts);

    opts->nlp_opts = ocp_nlp_opts_assign(config, dims, c_ptr);
    c_ptr += ocp_nlp_opts_calculate_size(config, dims);

    assert((char *) raw_memory + ocp_nlp_ddp_opts_calculate_size(config, dims) >= c_ptr);

    return opts;
}



void ocp_nlp_ddp_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_ddp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;

    // this first !!!
    ocp_nlp_opts_initialize_default(config, dims, nlp_opts);

    // DDP opts
    opts->max_iter = 20;
    opts->tol_stat = 1e-8;
    opts->tol_eq   = 1e-8;
    opts->tol_ineq = 1e-8;
    opts->tol_comp = 1e-8;
    opts->tol_zero_res = 1e-12;

    opts->ext_qp_res = 0;

    opts->qp_warm_start = 0;
    opts->warm_start_first_qp = false;
    opts->rti_phase = 0;
    opts->eval_residual_at_max_iter = false;

    opts->linesearch_eta = 1e-6;
    opts->linesearch_minimum_step_size = 1e-17;
    opts->linesearch_step_size_reduction_factor = 0.5;

    // overwrite default submodules opts

    // qp tolerance
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_stat", &opts->tol_stat);
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_eq", &opts->tol_eq);
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_ineq", &opts->tol_ineq);
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_comp", &opts->tol_comp);

    return;
}



void ocp_nlp_ddp_opts_update(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_ddp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_nlp_opts_update(config, dims, nlp_opts);

    return;
}



void ocp_nlp_ddp_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_ddp_opts *opts = (ocp_nlp_ddp_opts *) opts_;
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
        else if (!strcmp(field, "ext_qp_res"))
        {
            int* ext_qp_res = (int *) value;
            opts->ext_qp_res = *ext_qp_res;
        }
        else if (!strcmp(field, "warm_start_first_qp"))
        {
            bool* warm_start_first_qp = (bool *) value;
            opts->warm_start_first_qp = *warm_start_first_qp;
        }
        else if (!strcmp(field, "rti_phase"))
        {
            int* rti_phase = (int *) value;
            if (*rti_phase < 0 || *rti_phase > 0)
            {
                printf("\nerror: ocp_nlp_ddp_opts_set: invalid value for rti_phase field.");
                printf("possible values are: 0\n");
                exit(1);
            }
            opts->rti_phase = *rti_phase;
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



void ocp_nlp_ddp_opts_set_at_stage(void *config_, void *opts_, size_t stage, const char *field, void* value)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_ddp_opts *opts = (ocp_nlp_ddp_opts *) opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_nlp_opts_set_at_stage(config, nlp_opts, stage, field, value);

    return;

}

/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_ddp_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_ddp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_ddp_memory);

    // nlp mem
    size += ocp_nlp_memory_calculate_size(config, dims, nlp_opts);

    // stat
    int stat_m = opts->max_iter+1;
    int stat_n = 7;
    if (opts->ext_qp_res)
        stat_n += 4;
    size += stat_n*stat_m*sizeof(double);

    // matrices
    int nu_max = 0;
    int nx_max = 0;
    for (int i = 0; i <= N; i++)
    {
        nu_max = nu_max > nu[i] ? nu_max : nu[i];
        nx_max = nx_max > nx[i] ? nx_max : nx[i];
    }

    size += 1*blasfeo_memsize_dmat(nu_max, nx_max); // K_mat
    size += nu_max * nx_max * sizeof(double); // tmp_nu_times_nx

    size += 3*8;  // align
    size += 64;

    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_nlp_ddp_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_ddp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    // ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    // ocp_nlp_dynamics_config **dynamics = config->dynamics;
    // ocp_nlp_cost_config **cost = config->cost;
    // ocp_nlp_constraints_config **constraints = config->constraints;

    char *c_ptr = (char *) raw_memory;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;

    // initial align
    align_char_to(8, &c_ptr);

    ocp_nlp_ddp_memory *mem = (ocp_nlp_ddp_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_ddp_memory);

    align_char_to(8, &c_ptr);

    // nlp mem
    mem->nlp_mem = ocp_nlp_memory_assign(config, dims, nlp_opts, c_ptr);
    c_ptr += ocp_nlp_memory_calculate_size(config, dims, nlp_opts);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // dzduxt
    int nu_max = 0;
    int nx_max = 0;
    for (int i = 0; i <= N; i++)
    {
        nu_max = nu_max > nu[i] ? nu_max : nu[i];
        nx_max = nx_max > nx[i] ? nx_max : nx[i];
    }
    assign_and_advance_blasfeo_dmat_mem(nu_max, nx_max, &mem->K_mat, &c_ptr);

    // stat
    mem->stat = (double *) c_ptr;
    mem->stat_m = opts->max_iter+1;
    mem->stat_n = 7;
    if (opts->ext_qp_res)
        mem->stat_n += 4;
    c_ptr += mem->stat_m*mem->stat_n*sizeof(double);

    // tmp_nu_times_nx
    mem->tmp_nu_times_nx = (double *) c_ptr;
    c_ptr += nu_max*nx_max*sizeof(double);

    mem->status = ACADOS_READY;

    align_char_to(8, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_ddp_memory_calculate_size(config, dims, opts) >= c_ptr);

    return mem;
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_nlp_ddp_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_ddp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    acados_size_t size = 0;

    // ddp
    size += sizeof(ocp_nlp_ddp_workspace);

    // nlp
    size += ocp_nlp_workspace_calculate_size(config, dims, nlp_opts);

    if (opts->ext_qp_res)
    {
        // qp res
        size += ocp_qp_res_calculate_size(dims->qp_solver->orig_dims);

        // qp res ws
        size += ocp_qp_res_workspace_calculate_size(dims->qp_solver->orig_dims);
    }

    return size;
}



static void ocp_nlp_ddp_cast_workspace(ocp_nlp_config *config, ocp_nlp_dims *dims,
         ocp_nlp_ddp_opts *opts, ocp_nlp_ddp_memory *mem, ocp_nlp_ddp_workspace *work)
{
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    // ddp
    char *c_ptr = (char *) work;
    c_ptr += sizeof(ocp_nlp_ddp_workspace);

    // nlp
    work->nlp_work = ocp_nlp_workspace_assign(config, dims, nlp_opts, nlp_mem, c_ptr);
    c_ptr += ocp_nlp_workspace_calculate_size(config, dims, nlp_opts);

    if (opts->ext_qp_res)
    {
        // qp res
        work->qp_res = ocp_qp_res_assign(dims->qp_solver->orig_dims, c_ptr);
        c_ptr += ocp_qp_res_calculate_size(dims->qp_solver->orig_dims);

        // qp res ws
        work->qp_res_ws = ocp_qp_res_workspace_assign(dims->qp_solver->orig_dims, c_ptr);
        c_ptr += ocp_qp_res_workspace_calculate_size(dims->qp_solver->orig_dims);
    }

    assert((char *) work + ocp_nlp_ddp_workspace_calculate_size(config, dims, opts) >= c_ptr);

    return;
}


/************************************************
 * Helper functions
 ************************************************/
static void ocp_nlp_ddp_reset_timers(ocp_nlp_ddp_memory *mem)
{
    mem->time_qp_sol = 0.0;
    mem->time_qp_solver_call = 0.0;
    mem->time_qp_xcond = 0.0;
    mem->time_lin = 0.0;
    mem->time_reg = 0.0;
    mem->time_glob = 0.0;
    mem->time_sim = 0.0;
    mem->time_sim_la = 0.0;
    mem->time_sim_ad = 0.0;
}

static void ocp_nlp_ddp_compute_trial_iterate(ocp_nlp_config *config, ocp_nlp_dims *dims,
            ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_ddp_memory *ddp_mem,
            ocp_nlp_workspace *work, double alpha)
{
    /* computes trial iterate in tmp_nlp_out */
    int N = dims->N;
    // int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;
    int *nz = dims->nz;

    ocp_nlp_memory *mem = ddp_mem->nlp_mem;

    struct blasfeo_dvec *tmp_vec;
    ocp_nlp_out *tmp_nlp_out = work->tmp_nlp_out;
    ocp_qp_xcond_solver_config *xcond_solver_config = config->qp_solver;

    // compute x_0
    int i = 0;
    blasfeo_daxpy(nx[i], alpha, mem->qp_out->ux + i, nu[i],
                out->ux + i, nu[i], tmp_nlp_out->ux + i, nu[i]);

    // compute u_i, x_{i+1}
    for (i = 0; i < N; i++)
    {
        /* step in primal variables */
        // compute u   // (if i < N?)
        /* u_i = \bar{u}_i + alpha * k_i + K_i * (x_i - \bar{x}_i) */
        // get K
        xcond_solver_config->solver_get(xcond_solver_config, mem->qp_in, mem->qp_out, opts->qp_solver_opts, mem->qp_solver_mem, "K", i, ddp_mem->tmp_nu_times_nx, nu[i], nx[i]);
        blasfeo_pack_dmat(nu[i], nx[i], ddp_mem->tmp_nu_times_nx, nu[i], &ddp_mem->K_mat, 0, 0);

        // get k = tmp_nv;
        xcond_solver_config->solver_get(xcond_solver_config, mem->qp_in, mem->qp_out, opts->qp_solver_opts, mem->qp_solver_mem, "k", i, work->tmp_nv_double, nu[i], 1);
        blasfeo_pack_dvec(nu[i], work->tmp_nv_double, 1, &work->tmp_nv, 0);

        // compute delta_u = alpha * k_i + K_i * (x_i - \bar{x}_i)
        // tmp_nv[nu:] = (x_i - \bar{x}_i)
        blasfeo_daxpby(nx[i], -1.0, out->ux+i, nu[i], 1.0, tmp_nlp_out->ux+i, nu[i], &work->tmp_nv, nu[i]);
        blasfeo_dgemv_n(nu[i], nx[i], 1.0, &ddp_mem->K_mat, 0, 0, &work->tmp_nv, nu[i], alpha, &work->tmp_nv, 0, &work->tmp_nv, 0);
        blasfeo_daxpby(nu[i], 1.0, out->ux+i, 0, 1.0, &work->tmp_nv, 0, tmp_nlp_out->ux+i, 0);

        // evalutate dynamics
        // x_{i+1} = f_dyn_i(x_i, u_i)
        config->dynamics[i]->memory_set_ux_ptr(tmp_nlp_out->ux+i, mem->dynamics[i]);
        config->dynamics[i]->compute_fun(config->dynamics[i], dims->dynamics[i],
            in->dynamics[i], opts->dynamics[i], mem->dynamics[i], work->dynamics[i]);
        config->dynamics[i]->memory_set_ux_ptr(out->ux+i, mem->dynamics[i]);

        // f_dyn_i(x_i, u_i) - x_{i+1}
        // NOTE/TODO: store function output in dynamics module instead?
        tmp_vec = config->dynamics[i]->memory_get_fun_ptr(mem->dynamics[i]);
        blasfeo_daxpby(nx[i+1], 1.0, tmp_vec, 0, 1.0, out->ux+i+1, nu[i+1], tmp_nlp_out->ux+i+1, nu[i+1]);
    }

    for (i = 0; i < N+1; i++)
    {
        // update dual variables
        if (opts->full_step_dual)
        {
            blasfeo_dveccp(2*ni[i], mem->qp_out->lam+i, 0, tmp_nlp_out->lam+i, 0);
            if (i < N)
            {
                blasfeo_dveccp(nx[i+1], mem->qp_out->pi+i, 0, tmp_nlp_out->pi+i, 0);
            }
        }
        else
        {
            // update duals with alpha step
            blasfeo_dvecsc(2*ni[i], 1.0-alpha, tmp_nlp_out->lam+i, 0);
            blasfeo_daxpy(2*ni[i], alpha, mem->qp_out->lam+i, 0, tmp_nlp_out->lam+i, 0, tmp_nlp_out->lam+i, 0);
            if (i < N)
            {
                blasfeo_dvecsc(nx[i+1], 1.0-alpha, tmp_nlp_out->pi+i, 0);
                blasfeo_daxpy(nx[i+1], alpha, mem->qp_out->pi+i, 0, tmp_nlp_out->pi+i, 0, tmp_nlp_out->pi+i, 0);
            }
        }

        // linear update of algebraic variables using state and input sensitivity
        if (i < N)
        {
            // tmp_nlp_out->z = mem->z_alg + alpha * dzdux * qp_out->ux
            blasfeo_dgemv_t(nu[i]+nx[i], nz[i], alpha, mem->dzduxt+i, 0, 0,
                    mem->qp_out->ux+i, 0, 1.0, mem->z_alg+i, 0, tmp_nlp_out->z+i, 0);
        }
    }
}

/************************************************
 * output functions
 ************************************************/
static void print_iteration_header(){
    printf("%6s | %11s | %10s | %10s | %10s | %10s | %10s | %10s | %10s\n",
    "iter.",
    "objective",
    "res_eq",
    "res_stat",
    "alpha",
    "step_norm",
    "LM_reg.",
    "qp_status",
    "qp_iter");
}

static void print_iteration(double obj,
                     int iter_count,
                     double infeas,
                     double stationarity,
                     double alpha,
                     double step_norm,
                     double reg_param,
                     int qp_status,
                     int qp_iter)
{
    if ((iter_count % 10 == 0)){
        print_iteration_header();
    }
    printf("%6i | %11.4e | %10.4e | %10.4e | %10.4e | %10.4e | %10.4e | %10i | %10i\n",
    iter_count,
    obj,
    infeas,
    stationarity,
    alpha,
    step_norm,
    reg_param,
    qp_status,
    qp_iter);
}

/************************************************
 * termination criterion
 ************************************************/
static bool check_termination(int ddp_iter, ocp_nlp_res *nlp_res, ocp_nlp_ddp_memory *mem, ocp_nlp_ddp_opts *opts)
{
    // check for nans
    if (isnan(nlp_res->inf_norm_res_stat) || isnan(nlp_res->inf_norm_res_eq) ||
            isnan(nlp_res->inf_norm_res_ineq))
    {
        mem->status = ACADOS_NAN_DETECTED;
        if (opts->nlp_opts->print_level > 0)
        {
            printf("Stopped: NaN detected in iterate.\n");
        }
        return true;
    }

    // We do not need to check for the complementarity condition and for the
    // inequalities since we have an unconstrainted OCP
    if (nlp_res->inf_norm_res_eq < opts->tol_eq)
    { // Check that iterate must be dynamically feasible
        if (nlp_res->inf_norm_res_stat < opts->tol_stat)
        {// Check Stationarity
            mem->status = ACADOS_SUCCESS;
            if (opts->nlp_opts->print_level > 0)
            {
                printf("Optimal Solution found! Converged to KKT point.\n");
            }
            return true;
        }

        // Check for zero-residual solution of a least-squares problem
        if (opts->nlp_opts->with_adaptive_levenberg_marquardt && (mem->nlp_mem->cost_value < opts->tol_zero_res))
        {
            mem->status = ACADOS_SUCCESS;
            if (opts->nlp_opts->print_level > 0)
            {
                printf("Optimal Solution found! Converged To Zero Residual Solution.\n");
            }
            return true;
        }
    }

    // Check for small step
    if ((ddp_iter > 0) && (mem->step_norm < opts->tol_eq))
    {
        if (opts->nlp_opts->print_level > 0)
        {
            if (nlp_res->inf_norm_res_eq < opts->tol_eq)
            {
                printf("Stopped: Converged To Feasible Point. Step size is < tol_eq.\n");
            }
            else
            {
                printf("Stopped: Converged To Infeasible Point. Step size is < tol_eq.\n");
            }
        }
        mem->status = ACADOS_MINSTEP;
        return true;
    }

    // Check for maximum iterations
    if (ddp_iter >= opts->max_iter)
    {
        mem->status = ACADOS_MAXITER;
        if (opts->nlp_opts->print_level > 0){
            printf("Stopped: Maximum Iterations Reached.\n");
        }
        return true;
    }

    return false;
}


/************************************************
 * functions
 ************************************************/

// MAIN OPTIMIZATION ROUTINE
int ocp_nlp_ddp(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_)
{
    acados_timer timer0, timer1;
    acados_tic(&timer0);

    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_ddp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_ddp_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_res *nlp_res = nlp_mem->nlp_res;

    ocp_nlp_ddp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    ocp_qp_in *qp_in = nlp_mem->qp_in;
    ocp_qp_out *qp_out = nlp_mem->qp_out;

    // zero timers
    double tmp_time;
    ocp_nlp_ddp_reset_timers(mem);

    int N = dims->N;
    int ii;
    int qp_status = 0;
    int qp_iter = 0;
    mem->alpha = 0.0;
    mem->step_norm = 0.0;

#if defined(ACADOS_WITH_OPENMP)
    // backup number of threads
    int num_threads_bkp = omp_get_num_threads();
    // set number of threads
    // approximate_qp_matrices is parallelized
    omp_set_num_threads(opts->nlp_opts->num_threads);
#endif

    ocp_nlp_initialize_submodules(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);

    /************************************************
     * main ddp loop
     ************************************************/
    int ddp_iter = 0;
    double reg_param_memory = 0.0;
    bool infeasible_initial_guess = true;
    // bool evaluate_cost = true;
    if (nlp_opts->print_level > 0)
    {
        printf("'with_adaptive_levenberg_marquardt' option is set to: %s\n", opts->nlp_opts->with_adaptive_levenberg_marquardt?"true":"false");
    }

    for (; ddp_iter <= opts->max_iter; ddp_iter++)
    {
        // We always evaluate the residuals until the last iteration
        // If the option "eval_residual_at_max_iter" is set, then we will also
        // evaluate the data after the last iteration was performed
        if (ddp_iter != opts->max_iter || opts->eval_residual_at_max_iter)
        {
            /* Prepare the QP data */
            // linearize NLP, update QP matrices, and add Levenberg-Marquardt term
            acados_tic(&timer1);
            ocp_nlp_approximate_qp_matrices(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
            if (nlp_opts->with_adaptive_levenberg_marquardt || nlp_opts->globalization != FIXED_STEP)
            {
                ocp_nlp_get_cost_value_from_submodules(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
            }
            ocp_nlp_add_levenberg_marquardt_term(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, mem->alpha, ddp_iter);

            mem->time_lin += acados_toc(&timer1);

            // get timings from integrator
            for (ii=0; ii<N; ii++)
            {
                config->dynamics[ii]->memory_get(config->dynamics[ii], dims->dynamics[ii], mem->nlp_mem->dynamics[ii], "time_sim", &tmp_time);
                mem->time_sim += tmp_time;
                config->dynamics[ii]->memory_get(config->dynamics[ii], dims->dynamics[ii], mem->nlp_mem->dynamics[ii], "time_sim_la", &tmp_time);
                mem->time_sim_la += tmp_time;
                config->dynamics[ii]->memory_get(config->dynamics[ii], dims->dynamics[ii], mem->nlp_mem->dynamics[ii], "time_sim_ad", &tmp_time);
                mem->time_sim_ad += tmp_time;
            }

            // update QP rhs for DDP (step prim var, abs dual var)
            // NOTE: The ddp version of approximate does not exist!
            ocp_nlp_approximate_qp_vectors_sqp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);

            // compute nlp residuals
            ocp_nlp_res_compute(dims, nlp_in, nlp_out, nlp_res, nlp_mem);
            ocp_nlp_res_get_inf_norm(nlp_res, &nlp_out->inf_norm_res);
        }

        // save statistics
        if ((ddp_iter < mem->stat_m) & (ddp_iter >= 0))
        {
            mem->stat[mem->stat_n*ddp_iter+0] = nlp_res->inf_norm_res_stat;
            mem->stat[mem->stat_n*ddp_iter+1] = nlp_res->inf_norm_res_eq;
            mem->stat[mem->stat_n*ddp_iter+2] = nlp_mem->cost_value;
            mem->stat[mem->stat_n*ddp_iter+3] = reg_param_memory;
        }

        // Check if initial guess was infeasible
        if ((infeasible_initial_guess == true) && (nlp_res->inf_norm_res_eq > opts->tol_eq))
        {
            if (nlp_opts->print_level > 0)
            {
                printf("Initial guess was infeasible!\n");
            }
        }
        else
        {
            infeasible_initial_guess = false;
        }

        // Output
        if (nlp_opts->print_level > 0)
        {
            print_iteration(nlp_mem->cost_value, ddp_iter, nlp_res->inf_norm_res_eq, nlp_res->inf_norm_res_stat, mem->alpha, mem->step_norm, reg_param_memory, qp_status, qp_iter);
        }
        reg_param_memory = nlp_opts->levenberg_marquardt;

        // Termination
        if (check_termination(ddp_iter, nlp_res, mem, opts))
        {
#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif
            mem->ddp_iter = ddp_iter;
            mem->time_tot = acados_toc(&timer0);
            return mem->status;
        }

        // regularize Hessian
        acados_tic(&timer1);
        config->regularize->regularize(config->regularize, dims->regularize,
                                               nlp_opts->regularize, nlp_mem->regularize_mem);
        mem->time_reg += acados_toc(&timer1);

        /* solve QP */
        // (typically) no warm start at first iteration
        if (ddp_iter == 0 && !opts->warm_start_first_qp)
        {
            int tmp_int = 0;
            config->qp_solver->opts_set(config->qp_solver, nlp_opts->qp_solver_opts,
                                         "warm_start", &tmp_int);
        }
        // Show input to QP
        // if (nlp_opts->print_level > ddp_iter + 1)
        if (nlp_opts->print_level > 1)
        {
            printf("\n\nDDP: ocp_qp_in at iteration %d\n", ddp_iter + 1);
            print_ocp_qp_in(qp_in);
        }

        // solve qp
        acados_tic(&timer1);
        qp_status = qp_solver->evaluate(qp_solver, dims->qp_solver, qp_in, qp_out,
                                        nlp_opts->qp_solver_opts, nlp_mem->qp_solver_mem, nlp_work->qp_work);
        mem->time_qp_sol += acados_toc(&timer1);

        qp_solver->memory_get(qp_solver, nlp_mem->qp_solver_mem, "time_qp_solver_call", &tmp_time);
        mem->time_qp_solver_call += tmp_time;
        qp_solver->memory_get(qp_solver, nlp_mem->qp_solver_mem, "time_qp_xcond", &tmp_time);
        mem->time_qp_xcond += tmp_time;

        // compute correct dual solution in case of Hessian regularization
        acados_tic(&timer1);
        config->regularize->correct_dual_sol(config->regularize, dims->regularize,
                                             nlp_opts->regularize, nlp_mem->regularize_mem);
        mem->time_reg += acados_toc(&timer1);

        // restore default warm start
        if (ddp_iter==0)
        {
            config->qp_solver->opts_set(config->qp_solver, nlp_opts->qp_solver_opts,
                                        "warm_start", &opts->qp_warm_start);
        }

        // if (nlp_opts->print_level > ddp_iter + 1)
        if (nlp_opts->print_level > 1)
        {
            printf("\n\nDDP: ocp_qp_out at iteration %d\n", ddp_iter + 1);
            print_ocp_qp_out(qp_out);
        }

        qp_info *qp_info_;
        ocp_qp_out_get(qp_out, "qp_info", &qp_info_);
        qp_iter = qp_info_->num_iter;

        // save statistics of last qp solver call
        if (ddp_iter+1 < mem->stat_m)
        {
            mem->stat[mem->stat_n*(ddp_iter+1)+4] = qp_status;
            mem->stat[mem->stat_n*(ddp_iter+1)+5] = qp_iter;
        }

        // compute external QP residuals (for debugging)
        if (opts->ext_qp_res)
        {
            ocp_qp_res_compute(qp_in, qp_out, work->qp_res, work->qp_res_ws);
            if (ddp_iter+1 < mem->stat_m)
                ocp_qp_res_compute_nrm_inf(work->qp_res, mem->stat+(mem->stat_n*(ddp_iter+1)+7));
        }

        // exit conditions on QP status
        if ((qp_status!=ACADOS_SUCCESS) & (qp_status!=ACADOS_MAXITER))
        {
            // increment ddp_iter to return full statistics and improve output below.
            ddp_iter++;

#ifndef ACADOS_SILENT
            printf("\nQP solver returned error status %d in DDP iteration %d, QP iteration %d.\n",
                   qp_status, ddp_iter, qp_iter);
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

            mem->status = ACADOS_QP_FAILURE;
            mem->ddp_iter = ddp_iter;
            mem->time_tot = acados_toc(&timer0);

            return mem->status;
        }

        // Compute the optimal QP objective function value
        nlp_mem->qp_cost_value = ocp_nlp_ddp_compute_qp_objective_value(dims, qp_in, qp_out,nlp_work, nlp_mem);

        // Calculate step norm
        // res_comp
        mem->step_norm = 0.0;
        double tmp_norm = 0.0;
        for (int i = 0; i <= N; i++)
        {
            int sum = dims->nx[i]+dims->nu[i];
            blasfeo_dvecnrm_inf(sum, &qp_out->ux[i], 0, &tmp_norm);
            mem->step_norm = tmp_norm > mem->step_norm ? tmp_norm : mem->step_norm;
        }
        /* end solve QP */

        /* globalization */
        if (infeasible_initial_guess)
        {
            // Accept the forward simulation to get feasible initial guess
            mem->alpha = 1.0;  // full step to obtain feasible initial gues
            ocp_nlp_ddp_compute_trial_iterate(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work, mem->alpha);
            copy_ocp_nlp_out(dims, work->nlp_work->tmp_nlp_out, nlp_out);
            infeasible_initial_guess = false;
        }
        else
        {
            int linesearch_success = 1;
            // Do the globalization here: Either fixed step or Armijo line search
            acados_tic(&timer1);
            // NOTE on timings: currently all within globalization is accounted for within time_glob.
            //   QP solver times could be also attributed there alternatively. Cleanest would be to save them seperately.
            if (nlp_opts->globalization == FIXED_STEP)
            {
                // Set the given step length
                mem->alpha = nlp_opts->step_length;

                // update variables
                ocp_nlp_ddp_compute_trial_iterate(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work, mem->alpha);
            }
            else if (nlp_opts->globalization == MERIT_BACKTRACKING)
            {
                // do backtracking line search on objective function
                linesearch_success = ocp_nlp_ddp_backtracking_line_search(config, dims, nlp_in, nlp_out, mem, work, opts);
                // evaluate_cost = false; // since the cost was already evaluated in the line search
            }

            mem->stat[mem->stat_n*(ddp_iter+1)+6] = mem->alpha;
            // Copy new iterate to nlp_out
            if (linesearch_success == 1)
            {
                // in case line search fails, we do not want to copy trial iterates!
                copy_ocp_nlp_out(dims, work->nlp_work->tmp_nlp_out, nlp_out);
            }
            mem->time_glob += acados_toc(&timer1);
        }
    }  // end DDP loop

    if (nlp_opts->print_level > 0)
    {
        printf("Warning: The solver should never reach this part of the function!\n");
    }
    return mem->status;
}

double ocp_nlp_ddp_compute_qp_objective_value(ocp_nlp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out,
                ocp_nlp_workspace *nlp_work, ocp_nlp_memory *nlp_mem){

    // Compute the QP objective function value
    double qp_cost = 0.0;
    int i, nux;
    int N = dims->N;
    // Sum over stages 0 to N
    for (i = 0; i <= N; i++)
    {
        nux = dims->nx[i] + dims->nu[i];
        // Calculate 0.5* d.T H d
        blasfeo_dsymv_l(nux, 0.5, &qp_in->RSQrq[i], 0, 0, &qp_out->ux[i], 0, 0.0, &qp_out->ux[i], 0, &nlp_work->tmp_nlp_out->ux[i], 0);
        qp_cost += blasfeo_ddot(nux, &qp_out->ux[i], 0, &nlp_work->tmp_nlp_out->ux[i], 0);
        // Calculate g.T d
        qp_cost += blasfeo_ddot(nux, &qp_out->ux[i], 0, &qp_in->rqz[i], 0);
    }
    return qp_cost;
}


int ocp_nlp_ddp_backtracking_line_search(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *mem_, void *work_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    ocp_nlp_ddp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_nlp_ddp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    ocp_nlp_ddp_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;

    // evaluate the objective of the QP (as predicted reduction)
    // double qp_cost = compute_qp_cost
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
        ocp_nlp_ddp_compute_trial_iterate(config, dims, nlp_in, nlp_out, nlp_opts, mem, nlp_work, alpha);

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
        if (negative_ared <= fmin(-opts->linesearch_eta*alpha* fmax(pred, 0) + 1e-18, 0))
        {
            // IF step accepted: update x
            // reset evaluation point to SQP iterate
            // copy_ocp_nlp_out(dims, work->nlp_work->tmp_nlp_out, nlp_out);
            mem->alpha = alpha;
            nlp_mem->cost_value = trial_cost;
            return 1;
        }
        else
        {
            // Reduce step size
            alpha *= opts->linesearch_step_size_reduction_factor;
        }

        if (alpha < opts->linesearch_minimum_step_size)
        {
            printf("Linesearch: Step size gets too small. Increasing regularization.\n");
            mem->alpha = 0.0; // set to zero such that regularization is increased
            return 0;
        }
    }
}

void ocp_nlp_ddp_memory_reset_qp_solver(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
    void *opts_, void *mem_, void *work_)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_ddp_opts *opts = opts_;
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_ddp_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_ddp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    // printf("in ocp_nlp_ddp_memory_reset_qp_solver\n\n");
    config->qp_solver->memory_reset(qp_solver, dims->qp_solver,
        nlp_mem->qp_in, nlp_mem->qp_out, opts->nlp_opts->qp_solver_opts,
        nlp_mem->qp_solver_mem, nlp_work->qp_work);
}


int ocp_nlp_ddp_precompute(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_ddp_opts *opts = opts_;
    ocp_nlp_ddp_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    nlp_mem->workspace_size = ocp_nlp_workspace_calculate_size(config, dims, opts->nlp_opts);

    ocp_nlp_ddp_workspace *work = work_;
    ocp_nlp_ddp_cast_workspace(config, dims, opts, mem, work);
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    // sanity checks
    int N = dims->N;
    int nbx;
    int i = 0;
    config->constraints[i]->dims_get(config->constraints[i], dims->constraints[i],
                    "nbx", &nbx);
    if (nbx > dims->ni[0])
    {
        printf("ocp_nlp_precompute: nbx_0 > ni[0] not supported, got %d > %d\n", nbx, dims->ni[0]);
        exit(1);
    }

    for (i = 1; i < N+1; i++)
    {
        if (dims->ni[i] > 0)
        {
            printf("ocp_nlp_ddp: only initial state constraints are supported.");
            printf("Got ni[%d] = %d.\n", i, dims->ni[i]);
            exit(1);
        }
    }

    return ocp_nlp_precompute_common(config, dims, nlp_in, nlp_out, opts->nlp_opts, nlp_mem, nlp_work);
}


void ocp_nlp_ddp_eval_param_sens(void *config_, void *dims_, void *opts_, void *mem_, void *work_,
                                 char *field, int stage, int index, void *sens_nlp_out_)
{
    acados_timer timer0;
    acados_tic(&timer0);

    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_ddp_opts *opts = opts_;
    ocp_nlp_ddp_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_nlp_out *sens_nlp_out = sens_nlp_out_;

    ocp_nlp_ddp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    ocp_nlp_common_eval_param_sens(config, dims, opts->nlp_opts, nlp_mem, nlp_work,
                                 field, stage, index, sens_nlp_out);

    mem->time_solution_sensitivities = acados_toc(&timer0);

    return;
}


void ocp_nlp_ddp_eval_lagr_grad_p(void *config_, void *dims_, void *nlp_in_, void *opts_, void *mem_, void *work_,
                                 const char *field, void *lagr_grad_wrt_params)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_ddp_opts *opts = opts_;
    ocp_nlp_ddp_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    ocp_nlp_in *nlp_in = nlp_in_;

    ocp_nlp_ddp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    ocp_nlp_common_eval_lagr_grad_p(config, dims, nlp_in, opts->nlp_opts, nlp_mem, nlp_work,
                                 field, lagr_grad_wrt_params);

    return;
}


void ocp_nlp_ddp_get(void *config_, void *dims_, void *mem_, const char *field, void *return_value_)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_ddp_memory *mem = mem_;

    if (!strcmp("ddp_iter", field) || !strcmp("nlp_iter", field))
    {
        int *value = return_value_;
        *value = mem->ddp_iter;
    }
    else if (!strcmp("status", field))
    {
        int *value = return_value_;
        *value = mem->status;
    }
    else if (!strcmp("time_tot", field) || !strcmp("tot_time", field))
    {
        double *value = return_value_;
        *value = mem->time_tot;
    }
    else if (!strcmp("time_qp_sol", field) || !strcmp("time_qp", field))
    {
        double *value = return_value_;
        *value = mem->time_qp_sol;
    }
    else if (!strcmp("time_qp_solver", field) || !strcmp("time_qp_solver_call", field))
    {
        double *value = return_value_;
        *value = mem->time_qp_solver_call;
    }
    else if (!strcmp("time_qp_xcond", field))
    {
        double *value = return_value_;
        *value = mem->time_qp_xcond;
    }
    else if (!strcmp("time_lin", field))
    {
        double *value = return_value_;
        *value = mem->time_lin;
    }
    else if (!strcmp("time_reg", field))
    {
        double *value = return_value_;
        *value = mem->time_reg;
    }
    else if (!strcmp("time_glob", field))
    {
        double *value = return_value_;
        *value = mem->time_glob;
    }
    else if (!strcmp("time_solution_sensitivities", field))
    {
        double *value = return_value_;
        *value = mem->time_solution_sensitivities;
    }
    else if (!strcmp("time_sim", field))
    {
        double *value = return_value_;
        *value = mem->time_sim;
    }
    else if (!strcmp("time_sim_la", field))
    {
        double *value = return_value_;
        *value = mem->time_sim_la;
    }
    else if (!strcmp("time_sim_ad", field))
    {
        double *value = return_value_;
        *value = mem->time_sim_ad;
    }
    else if (!strcmp("time_preparation", field))
    {
        double *value = return_value_;
        *value = 0.0;
    }
    else if (!strcmp("time_feedback", field))
    {
        double *value = return_value_;
        *value = mem->time_tot;
    }
    else if (!strcmp("stat", field))
    {
        double **value = return_value_;
        *value = mem->stat;
    }
    else if (!strcmp("statistics", field))
    {
        int n_row = mem->stat_m<mem->ddp_iter+1 ? mem->stat_m : mem->ddp_iter+1;
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
    else if (!strcmp("nlp_mem", field))
    {
        void **value = return_value_;
        *value = mem->nlp_mem;
    }
    else if (!strcmp("qp_xcond_dims", field))
    {
        void **value = return_value_;
        *value = dims->qp_solver->xcond_dims;
    }
    else if (!strcmp("nlp_res", field))
    {
        ocp_nlp_res **value = return_value_;
        *value = mem->nlp_mem->nlp_res;
    }
    else if (!strcmp("qp_xcond_in", field))
    {
        void **value = return_value_;
        *value = mem->nlp_mem->qp_solver_mem->xcond_qp_in;
    }
    else if (!strcmp("qp_xcond_out", field))
    {
        void **value = return_value_;
        *value = mem->nlp_mem->qp_solver_mem->xcond_qp_out;
    }
    else if (!strcmp("qp_in", field))
    {
        void **value = return_value_;
        *value = mem->nlp_mem->qp_in;
    }
    else if (!strcmp("qp_out", field))
    {
        void **value = return_value_;
        *value = mem->nlp_mem->qp_out;
    }
    else if (!strcmp("qp_iter", field))
    {
        config->qp_solver->memory_get(config->qp_solver,
            mem->nlp_mem->qp_solver_mem, "iter", return_value_);
    }
    else if (!strcmp("qp_status", field))
    {
        config->qp_solver->memory_get(config->qp_solver,
            mem->nlp_mem->qp_solver_mem, "status", return_value_);
    }
    else if (!strcmp("res_stat", field))
    {
        double *value = return_value_;
        *value = mem->nlp_mem->nlp_res->inf_norm_res_stat;
    }
    else if (!strcmp("res_eq", field))
    {
        double *value = return_value_;
        *value = mem->nlp_mem->nlp_res->inf_norm_res_eq;
    }
    else if (!strcmp("res_ineq", field))
    {
        double *value = return_value_;
        *value = mem->nlp_mem->nlp_res->inf_norm_res_ineq;
    }
    else if (!strcmp("res_comp", field))
    {
        double *value = return_value_;
        *value = mem->nlp_mem->nlp_res->inf_norm_res_comp;
    }
    else if (!strcmp("cost_value", field))
    {
        double *value = return_value_;
        *value = mem->nlp_mem->cost_value;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_ddp_get\n", field);
        exit(1);
    }
}



void ocp_nlp_ddp_opts_get(void *config_, void *dims_, void *opts_,
                          const char *field, void *return_value_)
{
    ocp_nlp_ddp_opts *opts = opts_;

    if (!strcmp("nlp_opts", field))
    {
        void **value = return_value_;
        *value = opts->nlp_opts;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_ddp_opts_get\n", field);
        exit(1);
    }
}


void ocp_nlp_ddp_work_get(void *config_, void *dims_, void *work_,
                          const char *field, void *return_value_)
{
    ocp_nlp_ddp_workspace *work = work_;

    if (!strcmp("nlp_work", field))
    {
        void **value = return_value_;
        *value = work->nlp_work;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_ddp_work_get\n", field);
        exit(1);
    }
}


void ocp_nlp_ddp_terminate(void *config_, void *mem_, void *work_)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_ddp_memory *mem = mem_;
    ocp_nlp_ddp_workspace *work = work_;

    config->qp_solver->terminate(config->qp_solver, mem->nlp_mem->qp_solver_mem, work->nlp_work->qp_work);
}


void ocp_nlp_ddp_config_initialize_default(void *config_)
{
    ocp_nlp_config *config = (ocp_nlp_config *) config_;

    config->opts_calculate_size = &ocp_nlp_ddp_opts_calculate_size;
    config->opts_assign = &ocp_nlp_ddp_opts_assign;
    config->opts_initialize_default = &ocp_nlp_ddp_opts_initialize_default;
    config->opts_update = &ocp_nlp_ddp_opts_update;
    config->opts_set = &ocp_nlp_ddp_opts_set;
    config->opts_set_at_stage = &ocp_nlp_ddp_opts_set_at_stage;
    config->memory_calculate_size = &ocp_nlp_ddp_memory_calculate_size;
    config->memory_assign = &ocp_nlp_ddp_memory_assign;
    config->workspace_calculate_size = &ocp_nlp_ddp_workspace_calculate_size;
    config->evaluate = &ocp_nlp_ddp;
    config->memory_reset_qp_solver = &ocp_nlp_ddp_memory_reset_qp_solver;
    config->eval_param_sens = &ocp_nlp_ddp_eval_param_sens;
    config->eval_lagr_grad_p = &ocp_nlp_ddp_eval_lagr_grad_p;
    config->config_initialize_default = &ocp_nlp_ddp_config_initialize_default;
    config->precompute = &ocp_nlp_ddp_precompute;
    config->get = &ocp_nlp_ddp_get;
    config->opts_get = &ocp_nlp_ddp_opts_get;
    config->work_get = &ocp_nlp_ddp_work_get;
    config->terminate = &ocp_nlp_ddp_terminate;

    return;
}

