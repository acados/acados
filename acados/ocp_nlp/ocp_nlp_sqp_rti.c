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


#include "acados/ocp_nlp/ocp_nlp_sqp_rti.h"

// external
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
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

acados_size_t ocp_nlp_sqp_rti_opts_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_sqp_rti_opts);

    size += ocp_nlp_opts_calculate_size(config, dims);

    return size;
}



void *ocp_nlp_sqp_rti_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_sqp_rti_opts *opts = (ocp_nlp_sqp_rti_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_sqp_rti_opts);

    opts->nlp_opts = ocp_nlp_opts_assign(config, dims, c_ptr);
    c_ptr += ocp_nlp_opts_calculate_size(config, dims);

    assert((char *) raw_memory + ocp_nlp_sqp_rti_opts_calculate_size(config,
        dims) >= c_ptr);

    return opts;
}



void ocp_nlp_sqp_rti_opts_initialize_default(void *config_,
    void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    // ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    // ocp_nlp_dynamics_config **dynamics = config->dynamics;
    // ocp_nlp_constraints_config **constraints = config->constraints;

    // int ii;

    // int N = dims->N;

    // this first !!!
    ocp_nlp_opts_initialize_default(config, dims, nlp_opts);

    // SQP RTI opts
    opts->ext_qp_res = 0;
    opts->warm_start_first_qp = false;
    opts->rti_phase = 0;

    // overwrite default submodules opts

    // do not compute adjoint in dynamics and constraints
    // int compute_adj = 0;

    // // dynamics
    // for (ii = 0; ii < N; ii++)
    // {
    //     dynamics[ii]->opts_set(dynamics[ii],
    //         opts->nlp_opts->dynamics[ii], "compute_adj", &compute_adj);
    // }

    // // constraints
    // for (ii = 0; ii <= N; ii++)
    // {
    //     constraints[ii]->opts_set(constraints[ii],
    //         opts->nlp_opts->constraints[ii], "compute_adj", &compute_adj);
    // }

    return;
}



void ocp_nlp_sqp_rti_opts_update(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_nlp_opts_update(config, dims, nlp_opts);

    return;
}



void ocp_nlp_sqp_rti_opts_set(void *config_, void *opts_,
    const char *field, void* value)
{
    ocp_nlp_sqp_rti_opts *opts = (ocp_nlp_sqp_rti_opts *) opts_;
    ocp_nlp_config *config = config_;
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
        if (!strcmp(field, "ext_qp_res"))
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
            if (*rti_phase < 0 || *rti_phase > 2) {
                printf("\nerror: ocp_nlp_sqp_opts_set: invalid value for rti_phase field."); 
                printf("possible values are: 0, 1, 2\n");
                exit(1);
            } else opts->rti_phase = *rti_phase;
        }
        else
        {
            ocp_nlp_opts_set(config, nlp_opts, field, value);
        }
    }

    return;

}



void ocp_nlp_sqp_rti_opts_set_at_stage(void *config_, void *opts_, size_t stage, const char *field, void* value)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = (ocp_nlp_sqp_rti_opts *) opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_nlp_opts_set_at_stage(config, nlp_opts, stage, field, value);
}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_sqp_rti_memory_calculate_size(void *config_,
    void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    // ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    // ocp_nlp_dynamics_config **dynamics = config->dynamics;
    // ocp_nlp_cost_config **cost = config->cost;
    // ocp_nlp_constraints_config **constraints = config->constraints;

    // int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;
    // int *nz = dims->nz;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_sqp_rti_memory);

    // nlp mem
    size += ocp_nlp_memory_calculate_size(config, dims, nlp_opts);

    // stat
    int stat_m = 1+1;
    int stat_n = 2;
    if (opts->ext_qp_res)
        stat_n += 4;
    size += stat_n*stat_m*sizeof(double);

    size += 8;  // initial align

    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_nlp_sqp_rti_memory_assign(void *config_, void *dims_,
    void *opts_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    // ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    // ocp_nlp_dynamics_config **dynamics = config->dynamics;
    // ocp_nlp_cost_config **cost = config->cost;
    // ocp_nlp_constraints_config **constraints = config->constraints;

    char *c_ptr = (char *) raw_memory;

    // int ii;

    // int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;
    // int *nz = dims->nz;

    // initial align
    align_char_to(8, &c_ptr);

    ocp_nlp_sqp_rti_memory *mem = (ocp_nlp_sqp_rti_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_sqp_rti_memory);

    // nlp mem
    mem->nlp_mem = ocp_nlp_memory_assign(config, dims, nlp_opts, c_ptr);
    c_ptr += ocp_nlp_memory_calculate_size(config, dims, nlp_opts);

    // stat
    mem->stat = (double *) c_ptr;
    mem->stat_m = 1+1;
    mem->stat_n = 2;
    if (opts->ext_qp_res)
        mem->stat_n += 4;
    c_ptr += mem->stat_m*mem->stat_n*sizeof(double);

    mem->status = ACADOS_READY;

    assert((char *) raw_memory+ocp_nlp_sqp_rti_memory_calculate_size(
        config, dims, opts) >= c_ptr);

    return mem;
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_nlp_sqp_rti_workspace_calculate_size(void *config_,
    void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    acados_size_t size = 0;

    // sqp
    size += sizeof(ocp_nlp_sqp_rti_workspace);

    // nlp
    size += ocp_nlp_workspace_calculate_size(config, dims, nlp_opts);

    // qp in
    size += ocp_qp_in_calculate_size(dims->qp_solver->orig_dims);

    // qp out
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



static void ocp_nlp_sqp_rti_cast_workspace(
    ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_sqp_rti_opts *opts, 
    ocp_nlp_sqp_rti_memory *mem, ocp_nlp_sqp_rti_workspace *work)
{
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    // sqp
    char *c_ptr = (char *) work;
    c_ptr += sizeof(ocp_nlp_sqp_rti_workspace);

    // nlp
    work->nlp_work = ocp_nlp_workspace_assign(
        config, dims, nlp_opts, nlp_mem, c_ptr);
    c_ptr += ocp_nlp_workspace_calculate_size(config, dims, nlp_opts);

    // qp in
    work->tmp_qp_in = ocp_qp_in_assign(dims->qp_solver->orig_dims, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(dims->qp_solver->orig_dims);

    // qp out
    work->tmp_qp_out = ocp_qp_out_assign(dims->qp_solver->orig_dims, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(dims->qp_solver->orig_dims);

    if (opts->ext_qp_res)
    {
        // qp res
        work->qp_res = ocp_qp_res_assign(dims->qp_solver->orig_dims, c_ptr);
        c_ptr += ocp_qp_res_calculate_size(dims->qp_solver->orig_dims);

        // qp res ws
        work->qp_res_ws = ocp_qp_res_workspace_assign(
            dims->qp_solver->orig_dims, c_ptr);
        c_ptr += ocp_qp_res_workspace_calculate_size(
            dims->qp_solver->orig_dims);
    }

    assert((char *) work + ocp_nlp_sqp_rti_workspace_calculate_size(config,
        dims, opts) >= c_ptr);

    return;
}



/************************************************
 * functions
 ************************************************/

int ocp_nlp_sqp_rti(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
    void *opts_, void *mem_, void *work_)
{
    ocp_nlp_sqp_rti_memory *mem = mem_;
    
    // zero timers
    acados_timer timer0;
    double total_time = 0.0;
    mem->time_tot = 0.0;

    ocp_nlp_sqp_rti_opts *nlp_opts = opts_;
    int rti_phase = nlp_opts->rti_phase; 

    acados_tic(&timer0);
    switch(rti_phase) 
    {
        
        // perform preparation and feedback rti_phase
        case 0:
            ocp_nlp_sqp_rti_preparation_step(
                config_, dims_, nlp_in_, nlp_out_, opts_, mem_, work_);

            ocp_nlp_sqp_rti_feedback_step(
                config_, dims_, nlp_in_, nlp_out_, opts_, mem_, work_);

            break;

        // perform preparation rti_phase
        case 1:
            ocp_nlp_sqp_rti_preparation_step(
                config_, dims_, nlp_in_, nlp_out_, opts_, mem_, work_);

            break;

        // perform feedback rti_phase
        case 2:
            ocp_nlp_sqp_rti_feedback_step(
                config_, dims_, nlp_in_, nlp_out_, opts_, mem_, work_);

            break;
    }

    total_time += acados_toc(&timer0);
    mem->time_tot = total_time;

    return mem->status;

}



void ocp_nlp_sqp_rti_preparation_step(void *config_, void *dims_,
    void *nlp_in_, void *nlp_out_, void *opts_, void *mem_, void *work_)
{
    acados_timer timer1;

    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_sqp_rti_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    // ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;

    ocp_nlp_sqp_rti_workspace *work = work_;
    ocp_nlp_sqp_rti_cast_workspace(config, dims, opts, mem, work);
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    mem->time_lin = 0.0;
    mem->time_reg = 0.0;

#if defined(ACADOS_WITH_OPENMP)
    // backup number of threads
    int num_threads_bkp = omp_get_num_threads();
    // set number of threads
    omp_set_num_threads(opts->nlp_opts->num_threads);
#endif

    ocp_nlp_alias_memory_to_submodules(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);

    // initialize QP
    ocp_nlp_initialize_submodules(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);

    /* SQP body */
    int sqp_iter = 0;
    nlp_mem->sqp_iter = &sqp_iter;

    // linearizate NLP and update QP matrices
    acados_tic(&timer1);
    ocp_nlp_approximate_qp_matrices(config, dims, nlp_in,
        nlp_out, nlp_opts, nlp_mem, nlp_work);

    mem->time_lin += acados_toc(&timer1);

#if defined(ACADOS_WITH_OPENMP)
    // restore number of threads
    omp_set_num_threads(num_threads_bkp);
#endif
	
	return;

}



void ocp_nlp_sqp_rti_feedback_step(void *config_, void *dims_,
    void *nlp_in_, void *nlp_out_, void *opts_, void *mem_, void *work_)
{
    acados_timer timer1;

    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_sqp_rti_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;

    ocp_nlp_sqp_rti_workspace *work = work_;
    ocp_nlp_sqp_rti_cast_workspace(config, dims, opts, mem, work);
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    int qp_iter = 0;
    int qp_status = 0;
    double tmp_time;
    mem->time_qp_sol = 0.0;
    mem->time_qp_solver_call = 0.0;
    mem->time_qp_xcond = 0.0;
    mem->time_glob = 0.0;

    // embed initial value (this actually updates all bounds at stage 0...)
    ocp_nlp_embed_initial_value(config, dims, nlp_in,
        nlp_out, nlp_opts, nlp_mem, nlp_work);

    // update QP rhs for SQP (step prim var, abs dual var)
    ocp_nlp_approximate_qp_vectors_sqp(config, dims, nlp_in,
        nlp_out, nlp_opts, nlp_mem, nlp_work);

    // regularize Hessian
    acados_tic(&timer1);
    config->regularize->regularize_hessian(config->regularize,
        dims->regularize, opts->nlp_opts->regularize, nlp_mem->regularize_mem);
    mem->time_reg += acados_toc(&timer1);

    if (nlp_opts->print_level > 0) {
        printf("\n------- qp_in --------\n");
        print_ocp_qp_in(nlp_mem->qp_in);
    }

    if (!opts->warm_start_first_qp)
    {
        int tmp_int = 0;
        config->qp_solver->opts_set(config->qp_solver,
            opts->nlp_opts->qp_solver_opts, "warm_start", &tmp_int);
    }

    // solve qp
    acados_tic(&timer1);
    qp_status = qp_solver->evaluate(qp_solver, dims->qp_solver,
        nlp_mem->qp_in, nlp_mem->qp_out, opts->nlp_opts->qp_solver_opts,
        nlp_mem->qp_solver_mem, nlp_work->qp_work);

    mem->time_qp_sol += acados_toc(&timer1);

    qp_solver->memory_get(qp_solver, nlp_mem->qp_solver_mem, "time_qp_solver_call", &tmp_time);
    mem->time_qp_solver_call += tmp_time;
    qp_solver->memory_get(qp_solver, nlp_mem->qp_solver_mem, "time_qp_xcond", &tmp_time);
    mem->time_qp_xcond += tmp_time;

    // compute correct dual solution in case of Hessian regularization
    acados_tic(&timer1);
    config->regularize->correct_dual_sol(config->regularize,
        dims->regularize, opts->nlp_opts->regularize, nlp_mem->regularize_mem);

    mem->time_reg += acados_toc(&timer1);

    // TODO move into QP solver memory ???
    qp_info *qp_info_;
    ocp_qp_out_get(nlp_mem->qp_out, "qp_info", &qp_info_);
    qp_iter = qp_info_->num_iter;

    // compute external QP residuals (for debugging)
    if (opts->ext_qp_res)
    {
        ocp_qp_res_compute(nlp_mem->qp_in, nlp_mem->qp_out,
            work->qp_res, work->qp_res_ws);

        ocp_qp_res_compute_nrm_inf(work->qp_res, mem->stat+(mem->stat_n*1+2));
//      printf("\nsqp_iter %d, res %e %e %e %e\n", sqp_iter,
//      inf_norm_qp_res[0], inf_norm_qp_res[1],
//      inf_norm_qp_res[2], inf_norm_qp_res[3]);
    }

    // printf("\n------- qp_out (sqp iter %d) ---------\n", sqp_iter);
    // print_ocp_qp_out(nlp_mem->qp_out);
    // exit(1);

    // save statistics
    mem->stat[mem->stat_n*1+0] = qp_status;
    mem->stat[mem->stat_n*1+1] = qp_iter;

    if ((qp_status!=ACADOS_SUCCESS) & (qp_status!=ACADOS_MAXITER))
    {
        //   print_ocp_qp_in(mem->qp_in);
#ifndef ACADOS_SILENT
        printf("\nSQP_RTI: QP solver returned error status %d QP iteration %d.\n",
                qp_status, qp_iter);
#endif
        if (nlp_opts->print_level > 0)
        {
            printf("\n Failed to solve the following QP:\n");
            print_ocp_qp_in(nlp_mem->qp_in);
        }
        mem->status = ACADOS_QP_FAILURE;
        return;
    }


    // globalization
    acados_tic(&timer1);
    double alpha = ocp_nlp_line_search(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, 0);
    mem->time_glob += acados_toc(&timer1);

    // update variables
    ocp_nlp_update_variables_sqp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, alpha);

    // ocp_nlp_dims_print(nlp_out->dims);
    // ocp_nlp_out_print(nlp_out);
    // exit(1);

    // print_ocp_qp_in(mem->qp_in);

    mem->status = ACADOS_SUCCESS;

}


void ocp_nlp_sqp_rti_memory_reset_qp_solver(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
    void *opts_, void *mem_, void *work_)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_sqp_rti_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_sqp_rti_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    // printf("in ocp_nlp_sqp_rti_memory_reset_qp_solver\n\n");
    config->qp_solver->memory_reset(qp_solver, dims->qp_solver,
        nlp_mem->qp_in, nlp_mem->qp_out, opts->nlp_opts->qp_solver_opts,
        nlp_mem->qp_solver_mem, nlp_work->qp_work);
}


int ocp_nlp_sqp_rti_precompute(void *config_, void *dims_, void *nlp_in_, 
    void *nlp_out_, void *opts_, void *mem_, void *work_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;
    ocp_nlp_sqp_rti_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    // ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    ocp_nlp_sqp_rti_workspace *work = work_;
    ocp_nlp_sqp_rti_cast_workspace(config, dims, opts, mem, work);
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    int N = dims->N;
    int status = ACADOS_SUCCESS;

    int ii;

    // TODO(giaf) flag to enable/disable checks
    for (ii = 0; ii <= N; ii++)
    {
        int module_val;
        config->constraints[ii]->dims_get(config->constraints[ii],
            dims->constraints[ii], "ns", &module_val);
        if (dims->ns[ii] != module_val)
        {
            printf("ocp_nlp_sqp_rti_precompute: inconsistent dimension ns \
                for stage %d with constraint module, got %d, module: %d.",
                ii, dims->ns[ii], module_val);
            exit(1);
        }
    }

    // precompute
    for (ii = 0; ii < N; ii++)
    {
        // set T
        config->dynamics[ii]->model_set(config->dynamics[ii],
            dims->dynamics[ii], nlp_in->dynamics[ii], "T", nlp_in->Ts+ii);

        // dynamics precompute
        status = config->dynamics[ii]->precompute(config->dynamics[ii],
            dims->dynamics[ii], nlp_in->dynamics[ii],
            opts->nlp_opts->dynamics[ii],
            nlp_mem->dynamics[ii],
            nlp_work->dynamics[ii]);

        if (status != ACADOS_SUCCESS)
            return status;
    }

    return status;
}



void ocp_nlp_sqp_rti_eval_param_sens(void *config_, void *dims_, void *opts_,
    void *mem_, void *work_, char *field, int stage, int index,
    void *sens_nlp_out_)
{
    acados_timer timer0;
    acados_tic(&timer0);

    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;
    ocp_nlp_sqp_rti_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_nlp_out *sens_nlp_out = sens_nlp_out_;

    ocp_nlp_sqp_rti_workspace *work = work_;
    ocp_nlp_sqp_rti_cast_workspace(config, dims, opts, mem, work);
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    d_ocp_qp_copy_all(nlp_mem->qp_in, work->tmp_qp_in);
    d_ocp_qp_set_rhs_zero(work->tmp_qp_in);

    double one = 1.0;

    if ((!strcmp("ex", field)) & (stage==0))
    {
        d_ocp_qp_set_el("lbx", stage, index, &one, work->tmp_qp_in);
        d_ocp_qp_set_el("ubx", stage, index, &one, work->tmp_qp_in);

//        d_ocp_qp_print(work->tmp_qp_in->dim, work->tmp_qp_in);

        config->qp_solver->eval_sens(config->qp_solver, dims->qp_solver,
            work->tmp_qp_in, work->tmp_qp_out, opts->nlp_opts->qp_solver_opts,
            nlp_mem->qp_solver_mem, nlp_work->qp_work);

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
            blasfeo_dveccp(nv[i], work->tmp_qp_out->ux + i, 0,
                sens_nlp_out->ux + i, 0);

            if (i < N)
                blasfeo_dveccp(nx[i + 1], work->tmp_qp_out->pi + i, 0,
                    sens_nlp_out->pi + i, 0);

            blasfeo_dveccp(2 * ni[i], work->tmp_qp_out->lam + i, 0,
                sens_nlp_out->lam + i, 0);

            blasfeo_dveccp(2 * ni[i], work->tmp_qp_out->t + i, 0,
                sens_nlp_out->t + i, 0);

        }

    }
    else
    {
        printf("\nerror: field %s at stage %d not available in \
            ocp_nlp_sqp_rti_eval_param_sens\n", field, stage);

        exit(1);
    }
    mem->time_solution_sensitivities = acados_toc(&timer0);

    return;
}



// TODO rename memory_get ???
void ocp_nlp_sqp_rti_get(void *config_, void *dims_, void *mem_,
    const char *field, void *return_value_)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_sqp_rti_memory *mem = mem_;

    if (!strcmp("sqp_iter", field))
    {
        int *value = return_value_;
        *value = 1;
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
    else if (!strcmp("time_sim", field) || !strcmp("time_sim_ad", field) || !strcmp("time_sim_la", field))
    {
        double tmp = 0.0;
        double *ptr = return_value_;
        int N = dims->N;
        int ii;
        *ptr = 0.0;
        for (ii=0; ii<N; ii++)
        {
            config->dynamics[ii]->memory_get(config->dynamics[ii], dims->dynamics[ii], mem->nlp_mem->dynamics[ii], field, &tmp);
            *ptr += tmp;
        }
    }
    else if (!strcmp("time_solution_sensitivities", field))
    {
        double *value = return_value_;
        *value = mem->time_solution_sensitivities;
    }
    else if (!strcmp("stat", field))
    {
        double **value = return_value_;
        *value = mem->stat;
    }
    else if (!strcmp("statistics", field))
    {
        int n_row = 2;
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
        printf("\nerror: field %s not available in ocp_nlp_sqp_rti_get\n", field);
        exit(1);
    }
}


void ocp_nlp_sqp_rti_opts_get(void *config_, void *dims_, void *opts_,
                          const char *field, void *return_value_)
{
    // ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;

    if (!strcmp("nlp_opts", field))
    {
        void **value = return_value_;
        *value = opts->nlp_opts;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_sqp_rti_opts_get\n", field);
        exit(1);
    }

}


void ocp_nlp_sqp_rti_work_get(void *config_, void *dims_, void *work_,
                          const char *field, void *return_value_)
{
    // ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_workspace *work = work_;

    if (!strcmp("nlp_work", field))
    {
        void **value = return_value_;
        *value = work->nlp_work;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_sqp_rti_work_get\n", field);
        exit(1);
    }
}

void ocp_nlp_sqp_rti_config_initialize_default(void *config_)
{
    ocp_nlp_config *config = (ocp_nlp_config *) config_;

    config->opts_calculate_size = &ocp_nlp_sqp_rti_opts_calculate_size;
    config->opts_assign = &ocp_nlp_sqp_rti_opts_assign;
    config->opts_initialize_default = &ocp_nlp_sqp_rti_opts_initialize_default;
    config->opts_update = &ocp_nlp_sqp_rti_opts_update;
    config->opts_set = &ocp_nlp_sqp_rti_opts_set;
    config->opts_set_at_stage = &ocp_nlp_sqp_rti_opts_set_at_stage;
    config->memory_calculate_size = &ocp_nlp_sqp_rti_memory_calculate_size;
    config->memory_assign = &ocp_nlp_sqp_rti_memory_assign;
    config->workspace_calculate_size = &ocp_nlp_sqp_rti_workspace_calculate_size;
    config->evaluate = &ocp_nlp_sqp_rti;
    config->memory_reset_qp_solver = &ocp_nlp_sqp_rti_memory_reset_qp_solver;
    config->eval_param_sens = &ocp_nlp_sqp_rti_eval_param_sens;
    config->config_initialize_default = &ocp_nlp_sqp_rti_config_initialize_default;
    config->precompute = &ocp_nlp_sqp_rti_precompute;
    config->get = &ocp_nlp_sqp_rti_get;
    config->opts_get = &ocp_nlp_sqp_rti_opts_get;
    config->work_get = &ocp_nlp_sqp_rti_work_get;

    return;
}

