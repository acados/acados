/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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


#include "acados/ocp_nlp/ocp_nlp_sqp.h"

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



/************************************************
 * options
 ************************************************/

int ocp_nlp_sqp_opts_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    int size = 0;

    size += sizeof(ocp_nlp_sqp_opts);

    size += ocp_nlp_opts_calculate_size(config, dims);

    return size;
}



void *ocp_nlp_sqp_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_sqp_opts *opts = (ocp_nlp_sqp_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_sqp_opts);

    opts->nlp_opts = ocp_nlp_opts_assign(config, dims, c_ptr);
    c_ptr += ocp_nlp_opts_calculate_size(config, dims);

    assert((char *) raw_memory + ocp_nlp_sqp_opts_calculate_size(config, dims) >= c_ptr);

    return opts;
}



void ocp_nlp_sqp_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;

    // int ii;

    // this first !!!
    ocp_nlp_opts_initialize_default(config, dims, nlp_opts);

    // SQP opts
    opts->max_iter = 20;
    opts->tol_stat = 1e-8;
    opts->tol_eq   = 1e-8;
    opts->tol_ineq = 1e-8;
    opts->tol_comp = 1e-8;

    opts->ext_qp_res = 0;

    opts->qp_warm_start = 0;
    opts->warm_start_first_qp = false;
    opts->rti_phase = 0;
    opts->print_level = 0;

    // overwrite default submodules opts

    // qp tolerance
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_stat", &opts->tol_stat);
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_eq", &opts->tol_eq);
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_ineq", &opts->tol_ineq);
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_comp", &opts->tol_comp);

    return;
}



void ocp_nlp_sqp_opts_update(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_nlp_opts_update(config, dims, nlp_opts);

    return;
}



void ocp_nlp_sqp_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = (ocp_nlp_sqp_opts *) opts_;
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
//        config->qp_solver->opts_set(config->qp_solver, opts->qp_solver_opts, field+module_length+1, value);
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
            if (*rti_phase < 0 || *rti_phase > 0) {
                printf("\nerror: ocp_nlp_sqp_opts_set: invalid value for rti_phase field."); 
                printf("possible values are: 0\n");
                exit(1);
            } else opts->rti_phase = *rti_phase;
        }
        else if (!strcmp(field, "print_level"))
        {
            int* print_level = (int *) value;
            if (*print_level < 0)
            {
                printf("\nerror: ocp_nlp_sqp_opts_set: invalid value for print_level field, need int >=0, got %d.", *print_level);
                exit(1);
            }
            opts->print_level = *print_level;
        }
        else
        {
            ocp_nlp_opts_set(config, nlp_opts, field, value);
//            printf("\nerror: ocp_nlp_sqp_opts_set: wrong field: %s\n", field);
//            exit(1);
        }
    }

    return;

}



void ocp_nlp_sqp_opts_set_at_stage(void *config_, void *opts_, int stage, const char *field, void* value)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = (ocp_nlp_sqp_opts *) opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_nlp_opts_set_at_stage(config, nlp_opts, stage, field, value);

    return;

}



/************************************************
 * memory
 ************************************************/

int ocp_nlp_sqp_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    // int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;
    // int *nz = dims->nz;

    int size = 0;

    size += sizeof(ocp_nlp_sqp_memory);

    // nlp res
    size += ocp_nlp_res_calculate_size(dims);

    // nlp mem
    size += ocp_nlp_memory_calculate_size(config, dims, nlp_opts);

    // stat
    int stat_m = opts->max_iter+1;
    int stat_n = 6;
    if (opts->ext_qp_res)
        stat_n += 4;
    size += stat_n*stat_m*sizeof(double);

    size += 8;  // initial align

    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_nlp_sqp_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    // ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    // ocp_nlp_dynamics_config **dynamics = config->dynamics;
    // ocp_nlp_cost_config **cost = config->cost;
    // ocp_nlp_constraints_config **constraints = config->constraints;

    char *c_ptr = (char *) raw_memory;

    // int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;
    // int *nz = dims->nz;

    // initial align
    align_char_to(8, &c_ptr);

    ocp_nlp_sqp_memory *mem = (ocp_nlp_sqp_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_sqp_memory);

    // nlp res
    mem->nlp_res = ocp_nlp_res_assign(dims, c_ptr);
    c_ptr += mem->nlp_res->memsize;

    // nlp mem
    mem->nlp_mem = ocp_nlp_memory_assign(config, dims, nlp_opts, c_ptr);
    c_ptr += ocp_nlp_memory_calculate_size(config, dims, nlp_opts);

    // stat
    mem->stat = (double *) c_ptr;
    mem->stat_m = opts->max_iter+1;
    mem->stat_n = 6;
    if (opts->ext_qp_res)
        mem->stat_n += 4;
    c_ptr += mem->stat_m*mem->stat_n*sizeof(double);

    mem->status = ACADOS_READY;

    assert((char *) raw_memory + ocp_nlp_sqp_memory_calculate_size(config, dims, opts) >= c_ptr);

    return mem;
}



/************************************************
 * workspace
 ************************************************/

int ocp_nlp_sqp_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    int size = 0;

    // sqp
    size += sizeof(ocp_nlp_sqp_workspace);

    // nlp
    size += ocp_nlp_workspace_calculate_size(config, dims, nlp_opts);

    // tmp qp in
    size += ocp_qp_in_calculate_size(dims->qp_solver->orig_dims);

    // tmp qp out
    size += ocp_qp_out_calculate_size(dims->qp_solver->orig_dims);

    if (opts->ext_qp_res)
    {
        // qp res
        size += ocp_qp_res_calculate_size(dims->qp_solver->orig_dims);

        // qp res ws
        size += ocp_qp_res_workspace_calculate_size(dims->qp_solver->orig_dims);
    }

    return size;
}



static void ocp_nlp_sqp_cast_workspace(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_sqp_opts *opts, ocp_nlp_sqp_memory *mem, ocp_nlp_sqp_workspace *work)
{
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    // sqp
    char *c_ptr = (char *) work;
    c_ptr += sizeof(ocp_nlp_sqp_workspace);

    // nlp
    work->nlp_work = ocp_nlp_workspace_assign(config, dims, nlp_opts, nlp_mem, c_ptr);
    c_ptr += ocp_nlp_workspace_calculate_size(config, dims, nlp_opts);

    // tmp qp in
    work->tmp_qp_in = ocp_qp_in_assign(dims->qp_solver->orig_dims, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(dims->qp_solver->orig_dims);

    // tmp qp out
    work->tmp_qp_out = ocp_qp_out_assign(dims->qp_solver->orig_dims, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(dims->qp_solver->orig_dims);

    if (opts->ext_qp_res)
    {
        // qp res
        work->qp_res = ocp_qp_res_assign(dims->qp_solver->orig_dims, c_ptr);
        c_ptr += ocp_qp_res_calculate_size(dims->qp_solver->orig_dims);

        // qp res ws
        work->qp_res_ws = ocp_qp_res_workspace_assign(dims->qp_solver->orig_dims, c_ptr);
        c_ptr += ocp_qp_res_workspace_calculate_size(dims->qp_solver->orig_dims);
    }

    assert((char *) work + ocp_nlp_sqp_workspace_calculate_size(config, dims, opts) >= c_ptr);

    return;
}



/************************************************
 * functions
 ************************************************/

int ocp_nlp_sqp(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_)
{

    acados_timer timer0, timer1;

    acados_tic(&timer0);

    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_sqp_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;

    ocp_nlp_sqp_workspace *work = work_;
    ocp_nlp_sqp_cast_workspace(config, dims, opts, mem, work);
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    // zero timers
    double total_time = 0.0;
	double tmp_time;
    mem->time_qp_sol = 0.0;
    mem->time_qp_solver_call = 0.0;
    mem->time_qp_xcond = 0.0;
    mem->time_lin = 0.0;
    mem->time_reg = 0.0;
    mem->time_tot = 0.0;

    int N = dims->N;

    int ii;

    int qp_iter = 0;
    int qp_status = 0;

#if defined(ACADOS_WITH_OPENMP)
    // backup number of threads
    int num_threads_bkp = omp_get_num_threads();
    // set number of threads
    omp_set_num_threads(opts->nlp_opts->num_threads);
    #pragma omp parallel
    { // beginning of parallel region
#endif

    // alias to dynamics_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for
#endif
    for (ii = 0; ii < N; ii++)
    {
        config->dynamics[ii]->memory_set_ux_ptr(nlp_out->ux+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_tmp_ux_ptr(nlp_work->tmp_nlp_out->ux+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_ux1_ptr(nlp_out->ux+ii+1, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_tmp_ux1_ptr(nlp_work->tmp_nlp_out->ux+ii+1, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_pi_ptr(nlp_out->pi+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_tmp_pi_ptr(nlp_work->tmp_nlp_out->pi+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_BAbt_ptr(nlp_mem->qp_in->BAbt+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_RSQrq_ptr(nlp_mem->qp_in->RSQrq+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_dzduxt_ptr(nlp_mem->dzduxt+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_sim_guess_ptr(nlp_mem->sim_guess+ii, nlp_mem->set_sim_guess+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_z_alg_ptr(nlp_mem->z_alg+ii, nlp_mem->dynamics[ii]);
    }

    // alias to cost_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for
#endif
    for (ii = 0; ii <= N; ii++)
    {
        config->cost[ii]->memory_set_ux_ptr(nlp_out->ux+ii, nlp_mem->cost[ii]);
        config->cost[ii]->memory_set_tmp_ux_ptr(nlp_work->tmp_nlp_out->ux+ii, nlp_mem->cost[ii]);
        config->cost[ii]->memory_set_z_alg_ptr(nlp_mem->z_alg+ii, nlp_mem->cost[ii]);
        config->cost[ii]->memory_set_dzdux_tran_ptr(nlp_mem->dzduxt+ii, nlp_mem->cost[ii]);
        config->cost[ii]->memory_set_RSQrq_ptr(nlp_mem->qp_in->RSQrq+ii, nlp_mem->cost[ii]);
        config->cost[ii]->memory_set_Z_ptr(nlp_mem->qp_in->Z+ii, nlp_mem->cost[ii]);
    }
    // alias to constraints_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for
#endif
    for (ii = 0; ii <= N; ii++)
    {
        config->constraints[ii]->memory_set_ux_ptr(nlp_out->ux+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_tmp_ux_ptr(nlp_work->tmp_nlp_out->ux+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_lam_ptr(nlp_out->lam+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_tmp_lam_ptr(nlp_work->tmp_nlp_out->lam+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_z_alg_ptr(nlp_mem->z_alg+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_dzdux_tran_ptr(nlp_mem->dzduxt+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_DCt_ptr(nlp_mem->qp_in->DCt+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_RSQrq_ptr(nlp_mem->qp_in->RSQrq+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_idxb_ptr(nlp_mem->qp_in->idxb[ii], nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_idxs_ptr(nlp_mem->qp_in->idxs[ii], nlp_mem->constraints[ii]);
    }

    // alias to regularize memory
    config->regularize->memory_set_RSQrq_ptr(dims->regularize, nlp_mem->qp_in->RSQrq, nlp_mem->regularize_mem);
    config->regularize->memory_set_rq_ptr(dims->regularize, nlp_mem->qp_in->rqz, nlp_mem->regularize_mem);
    config->regularize->memory_set_BAbt_ptr(dims->regularize, nlp_mem->qp_in->BAbt, nlp_mem->regularize_mem);
    config->regularize->memory_set_b_ptr(dims->regularize, nlp_mem->qp_in->b, nlp_mem->regularize_mem);
    config->regularize->memory_set_idxb_ptr(dims->regularize, nlp_mem->qp_in->idxb, nlp_mem->regularize_mem);
    config->regularize->memory_set_DCt_ptr(dims->regularize, nlp_mem->qp_in->DCt, nlp_mem->regularize_mem);
    config->regularize->memory_set_ux_ptr(dims->regularize, nlp_mem->qp_out->ux, nlp_mem->regularize_mem);
    config->regularize->memory_set_pi_ptr(dims->regularize, nlp_mem->qp_out->pi, nlp_mem->regularize_mem);
    config->regularize->memory_set_lam_ptr(dims->regularize, nlp_mem->qp_out->lam, nlp_mem->regularize_mem);

    // copy sampling times into dynamics model
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for
#endif

    // NOTE(oj): this will lead in an error for irk_gnsf, T must be set in precompute;
    //    -> remove here and make sure precompute is called everywhere (e.g. Python interface).
    for (ii = 0; ii < N; ii++)
    {
        config->dynamics[ii]->model_set(config->dynamics[ii], dims->dynamics[ii],
                                         nlp_in->dynamics[ii], "T", nlp_in->Ts+ii);
    }

#if defined(ACADOS_WITH_OPENMP)
    } // end of parallel region
#endif

    // initialize QP
    ocp_nlp_initialize_qp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);

    // main sqp loop
    int sqp_iter = 0;
	nlp_mem->sqp_iter = &sqp_iter;

    for (; sqp_iter < opts->max_iter; sqp_iter++)
    {
		
        if (opts->print_level > 0)
        {
            printf("\n------- sqp iter %d (max_iter %d) --------\n", 
                sqp_iter, opts->max_iter);
            if (opts->print_level > sqp_iter + 1)
                print_ocp_qp_in(nlp_mem->qp_in);
        }

        // linearizate NLP and update QP matrices
        acados_tic(&timer1);
        ocp_nlp_approximate_qp_matrices(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
        mem->time_lin += acados_toc(&timer1);

        // update QP rhs for SQP (step prim var, abs dual var)
        ocp_nlp_approximate_qp_vectors_sqp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);

        // compute nlp residuals
        ocp_nlp_res_compute(dims, nlp_in, nlp_out, mem->nlp_res, nlp_mem);

        nlp_out->inf_norm_res = mem->nlp_res->inf_norm_res_g;
        nlp_out->inf_norm_res = (mem->nlp_res->inf_norm_res_b > nlp_out->inf_norm_res) ?
                                    mem->nlp_res->inf_norm_res_b :
                                    nlp_out->inf_norm_res;
        nlp_out->inf_norm_res = (mem->nlp_res->inf_norm_res_d > nlp_out->inf_norm_res) ?
                                    mem->nlp_res->inf_norm_res_d :
                                    nlp_out->inf_norm_res;
        nlp_out->inf_norm_res = (mem->nlp_res->inf_norm_res_m > nlp_out->inf_norm_res) ?
                                    mem->nlp_res->inf_norm_res_m :
                                    nlp_out->inf_norm_res;

        // save statistics
        if (sqp_iter < mem->stat_m)
        {
            mem->stat[mem->stat_n*sqp_iter+0] = mem->nlp_res->inf_norm_res_g;
            mem->stat[mem->stat_n*sqp_iter+1] = mem->nlp_res->inf_norm_res_b;
            mem->stat[mem->stat_n*sqp_iter+2] = mem->nlp_res->inf_norm_res_d;
            mem->stat[mem->stat_n*sqp_iter+3] = mem->nlp_res->inf_norm_res_m;
        }

        // exit conditions on residuals
        if ((mem->nlp_res->inf_norm_res_g < opts->tol_stat) &
            (mem->nlp_res->inf_norm_res_b < opts->tol_eq) &
            (mem->nlp_res->inf_norm_res_d < opts->tol_ineq) &
            (mem->nlp_res->inf_norm_res_m < opts->tol_comp))
        {
            // save sqp iterations number
            mem->sqp_iter = sqp_iter;
            nlp_out->sqp_iter = sqp_iter;

            // stop timer
            total_time += acados_toc(&timer0);

            // save time
            nlp_out->total_time = total_time;
            mem->time_tot = total_time;

#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif
            mem->status = ACADOS_SUCCESS;
            return mem->status;
        }


        // regularize Hessian
        acados_tic(&timer1);
        config->regularize->regularize_hessian(config->regularize, dims->regularize,
                                               opts->nlp_opts->regularize, nlp_mem->regularize_mem);
        mem->time_reg += acados_toc(&timer1);

        // (typically) no warm start at first iteration
        if (sqp_iter == 0 && !opts->warm_start_first_qp)
        {
            int tmp_int = 0;
            config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts,
                                         "warm_start", &tmp_int);
        }

        // solve qp
        acados_tic(&timer1);
        qp_status = qp_solver->evaluate(qp_solver, dims->qp_solver, nlp_mem->qp_in, nlp_mem->qp_out,
                                        opts->nlp_opts->qp_solver_opts, nlp_mem->qp_solver_mem, nlp_work->qp_work);
        mem->time_qp_sol += acados_toc(&timer1);

		qp_solver->memory_get(qp_solver, nlp_mem->qp_solver_mem, "time_qp_solver_call", &tmp_time);
		mem->time_qp_solver_call += tmp_time;
		qp_solver->memory_get(qp_solver, nlp_mem->qp_solver_mem, "time_qp_xcond", &tmp_time);
		mem->time_qp_xcond += tmp_time;

        // compute correct dual solution in case of Hessian regularization
        acados_tic(&timer1);
        config->regularize->correct_dual_sol(config->regularize, dims->regularize,
                                             opts->nlp_opts->regularize, nlp_mem->regularize_mem);
        mem->time_reg += acados_toc(&timer1);

        // restore default warm start
        if (sqp_iter==0)
        {
            config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts,
                                        "warm_start", &opts->qp_warm_start);
        }

        // TODO move into QP solver memory ???
        qp_info *qp_info_;
        ocp_qp_out_get(nlp_mem->qp_out, "qp_info", &qp_info_);
        nlp_out->qp_iter = qp_info_->num_iter;
        // printf("\nqp_iter = %d, sqp_iter = %d, max_sqp_iter = %d\n", nlp_out->qp_iter, sqp_iter, opts->max_iter);
        qp_iter = qp_info_->num_iter;

        // save statistics of last qp solver call
        if (sqp_iter+1 < mem->stat_m)
        {
            mem->stat[mem->stat_n*(sqp_iter+1)+4] = qp_status;
            mem->stat[mem->stat_n*(sqp_iter+1)+5] = qp_iter;
        }

        // compute external QP residuals (for debugging)
        if (opts->ext_qp_res)
        {
            ocp_qp_res_compute(nlp_mem->qp_in, nlp_mem->qp_out, work->qp_res, work->qp_res_ws);
            if (sqp_iter+1 < mem->stat_m)
                ocp_qp_res_compute_nrm_inf(work->qp_res, mem->stat+(mem->stat_n*(sqp_iter+1)+6));
        }


        if ((qp_status!=ACADOS_SUCCESS) & (qp_status!=ACADOS_MAXITER))
        {
            // print_ocp_qp_in(nlp_mem->qp_in);

            // save sqp iterations number
            mem->sqp_iter = sqp_iter;
            nlp_out->sqp_iter = sqp_iter;

            // stop timer
            total_time += acados_toc(&timer0);

            // save time
            mem->time_tot = total_time;
            nlp_out->total_time = total_time;

            printf("QP solver returned error status %d in iteration %d\n", qp_status, sqp_iter);
#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif

            if (opts->print_level > 1)
            {
                printf("\n Failed to solve the following QP:\n");
                if (opts->print_level > sqp_iter + 1)
                    print_ocp_qp_in(nlp_mem->qp_in);
            }

            mem->status = ACADOS_QP_FAILURE;
            return mem->status;
        }

        ocp_nlp_update_variables_sqp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);

        // ocp_nlp_dims_print(nlp_out->dims);
        // ocp_nlp_out_print(nlp_out);
        // exit(1);

        // ??? @rien
        //        for (int_t i = 0; i < N; i++)
        //        {
        //   ocp_nlp_dynamics_opts *dynamics_opts = opts->dynamics[i];
        //            sim_opts *opts = dynamics_opts->sim_solver;
        //            if (opts->scheme == NULL)
        //                continue;
        //            opts->sens_adj = (opts->scheme->type != exact);
        //            if (nlp_in->freezeSens) {
        //                // freeze inexact sensitivities after first SQP iteration !!
        //                opts->scheme->freeze = true;
        //            }
        //        }
        if (opts->print_level > 0)
        {
            printf("Residuals: stat: %e, eq: %e, ineq: %e, comp: %e.\n", mem->nlp_res->inf_norm_res_g,
                    mem->nlp_res->inf_norm_res_b, mem->nlp_res->inf_norm_res_d, mem->nlp_res->inf_norm_res_m );
        }

    }

    // stop timer
    total_time += acados_toc(&timer0);

    // ocp_nlp_out_print(nlp_out);

    // save sqp iterations number
    mem->sqp_iter = sqp_iter;
    nlp_out->sqp_iter = sqp_iter;

    // save time
    mem->time_tot = total_time;
    nlp_out->total_time = total_time;

    // maximum number of iterations reached
#if defined(ACADOS_WITH_OPENMP)
    // restore number of threads
    omp_set_num_threads(num_threads_bkp);
#endif
    mem->status = ACADOS_MAXITER;
    printf("\n ocp_nlp_sqp: maximum iterations reached\n");

    if (opts->print_level > 0)
    {
        printf("Residuals: stat: %e, eq: %e, ineq: %e, comp: %e.\n", mem->nlp_res->inf_norm_res_g,
            mem->nlp_res->inf_norm_res_b, mem->nlp_res->inf_norm_res_d, mem->nlp_res->inf_norm_res_m );
    }

    return mem->status;
}



int ocp_nlp_sqp_precompute(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_sqp_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    // ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    ocp_nlp_sqp_workspace *work = work_;
    ocp_nlp_sqp_cast_workspace(config, dims, opts, mem, work);
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    int N = dims->N;
    int status = ACADOS_SUCCESS;

    int ii;

    // TODO(all) add flag to enable/disable checks
    for (ii = 0; ii <= N; ii++)
    {
        int module_val;
        config->constraints[ii]->dims_get(config->constraints[ii], dims->constraints[ii], "ns", &module_val);
        if (dims->ns[ii] != module_val)
        {
            printf("ocp_nlp_sqp_precompute: inconsistent dimension ns for stage %d with constraint module, got %d, module: %d.",
                   ii, dims->ns[ii], module_val);
            exit(1);
        }
    }

    // precompute
    for (ii = 0; ii < N; ii++)
    {
        // set T
        config->dynamics[ii]->model_set(config->dynamics[ii], dims->dynamics[ii],
                                        nlp_in->dynamics[ii], "T", nlp_in->Ts+ii);
        // dynamics precompute
        status = config->dynamics[ii]->precompute(config->dynamics[ii], dims->dynamics[ii],
                                                nlp_in->dynamics[ii], opts->nlp_opts->dynamics[ii],
                                                nlp_mem->dynamics[ii], nlp_work->dynamics[ii]);
        if (status != ACADOS_SUCCESS)
            return status;
    }
    return status;
}



void ocp_nlp_sqp_eval_param_sens(void *config_, void *dims_, void *opts_, void *mem_, void *work_,
                                 char *field, int stage, int index, void *sens_nlp_out_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_sqp_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_nlp_out *sens_nlp_out = sens_nlp_out_;

    ocp_nlp_sqp_workspace *work = work_;
    ocp_nlp_sqp_cast_workspace(config, dims, opts, mem, work);
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    d_ocp_qp_copy_all(nlp_mem->qp_in, work->tmp_qp_in);
    d_ocp_qp_set_rhs_zero(work->tmp_qp_in);

    double one = 1.0;

    if ((!strcmp("ex", field)) & (stage==0))
    {
        d_ocp_qp_set_el("lbx", stage, index, &one, work->tmp_qp_in);
        d_ocp_qp_set_el("ubx", stage, index, &one, work->tmp_qp_in);

//        d_ocp_qp_print(work->tmp_qp_in->dim, work->tmp_qp_in);

        config->qp_solver->eval_sens(config->qp_solver, dims->qp_solver, work->tmp_qp_in, work->tmp_qp_out,
                               opts->nlp_opts->qp_solver_opts, nlp_mem->qp_solver_mem, nlp_work->qp_work);

//        d_ocp_qp_sol_print(work->tmp_qp_out->dim, work->tmp_qp_out);
//        exit(1);
        
        /* copy tmp_qp_out into sens_nlp_out */

        int i;

        int N = dims->N;
        int *nv = dims->nv;
        int *nx = dims->nx;
        // int *nu = dims->nu;
        int *ni = dims->ni;
        // int *nz = dims->nz;

        for (i = 0; i <= N; i++)
        {
            blasfeo_dveccp(nv[i], work->tmp_qp_out->ux + i, 0, sens_nlp_out->ux + i, 0);

            if (i < N)
                blasfeo_dveccp(nx[i + 1], work->tmp_qp_out->pi + i, 0, sens_nlp_out->pi + i, 0);

            blasfeo_dveccp(2 * ni[i], work->tmp_qp_out->lam + i, 0, sens_nlp_out->lam + i, 0);

            blasfeo_dveccp(2 * ni[i], work->tmp_qp_out->t + i, 0, sens_nlp_out->t + i, 0);

        }

    }
    else
    {
        printf("\nerror: field %s at stage %d not available in ocp_nlp_sqp_eval_param_sens\n", field, stage);
        exit(1);
    }

    return;
}



// TODO rename memory_get ???
void ocp_nlp_sqp_get(void *config_, void *dims_, void *mem_, const char *field, void *return_value_)
{
    ocp_nlp_config *config = config_;
	ocp_nlp_dims *dims = dims_;
    ocp_nlp_sqp_memory *mem = mem_;

    if (!strcmp("sqp_iter", field))
    {
        int *value = return_value_;
        *value = mem->sqp_iter;
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
    else if (!strcmp("time_sim", field) || !strcmp("time_sim_ad", field) || !strcmp("time_sim_la", field))
    {
		double tmp = 0.0;
		double *ptr = return_value_;
		int N = dims->N;
		int ii;
		for (ii=0; ii<N; ii++)
		{
			config->dynamics[ii]->memory_get(config->dynamics[ii], dims->dynamics[ii], mem->nlp_mem->dynamics[ii], field, &tmp);
			*ptr += tmp;
		}
	}
    else if (!strcmp("nlp_res", field))
    {
        ocp_nlp_res **value = return_value_;
        *value = mem->nlp_res;
    }
    else if (!strcmp("stat", field))
    {
        double **value = return_value_;
        *value = mem->stat;
    }
    else if (!strcmp("statistics", field))
    {
        int n_row = mem->stat_m<mem->sqp_iter+1 ? mem->stat_m : mem->sqp_iter+1;
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
		config->qp_solver->memory_get(config->qp_solver, mem->nlp_mem->qp_solver_mem, "iter", return_value_);
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_sqp_get\n", field);
        exit(1);
    }

}



void ocp_nlp_sqp_config_initialize_default(void *config_)
{
    ocp_nlp_config *config = (ocp_nlp_config *) config_;

    config->opts_calculate_size = &ocp_nlp_sqp_opts_calculate_size;
    config->opts_assign = &ocp_nlp_sqp_opts_assign;
    config->opts_initialize_default = &ocp_nlp_sqp_opts_initialize_default;
    config->opts_update = &ocp_nlp_sqp_opts_update;
    config->opts_set = &ocp_nlp_sqp_opts_set;
    config->opts_set_at_stage = &ocp_nlp_sqp_opts_set_at_stage;
    config->memory_calculate_size = &ocp_nlp_sqp_memory_calculate_size;
    config->memory_assign = &ocp_nlp_sqp_memory_assign;
    config->workspace_calculate_size = &ocp_nlp_sqp_workspace_calculate_size;
    config->evaluate = &ocp_nlp_sqp;
    config->eval_param_sens = &ocp_nlp_sqp_eval_param_sens;
    config->config_initialize_default = &ocp_nlp_sqp_config_initialize_default;
    config->precompute = &ocp_nlp_sqp_precompute;
    config->get = &ocp_nlp_sqp_get;

    return;
}
