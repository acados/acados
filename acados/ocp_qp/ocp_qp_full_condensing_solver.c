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


// external
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_full_condensing.h" // TODO remove !!!
#include "acados/ocp_qp/ocp_qp_full_condensing_solver.h" // TODO remove !!!
#include "acados/ocp_qp/ocp_qp_xcond_solver.h"
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"



/************************************************
 * memory
 ************************************************/

#if 0
int ocp_qp_full_condensing_solver_memory_calculate_size(void *config_, ocp_qp_xcond_solver_dims *dims,
                                                        void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
	ocp_qp_xcond_config *xcond = config->xcond;

//    ocp_qp_full_condensing_solver_opts *opts = (ocp_qp_full_condensing_solver_opts *) opts_;
    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;

    void *xcond_dims;
	xcond->dims_get(xcond, dims->xcond_dims, "xcond_dims", &xcond_dims);

    int size = 0;
    size += sizeof(ocp_qp_full_condensing_solver_memory);

//    size += ocp_qp_full_condensing_memory_calculate_size(dims, opts->xcond_opts);
    size += ocp_qp_full_condensing_memory_calculate_size(dims->xcond_dims, opts->xcond_opts);

    size += qp_solver->memory_calculate_size(qp_solver, xcond_dims, opts->qp_solver_opts);

    size += dense_qp_in_calculate_size(xcond_dims);
    size += dense_qp_out_calculate_size(xcond_dims);

    return size;
}



void *ocp_qp_full_condensing_solver_memory_assign(void *config_, ocp_qp_xcond_solver_dims *dims, void *opts_,
                                                  void *raw_memory)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
	ocp_qp_xcond_config *xcond = config->xcond;

//    ocp_qp_full_condensing_solver_opts *opts = (ocp_qp_full_condensing_solver_opts *) opts_;
    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;

    void *xcond_dims;
	xcond->dims_get(xcond, dims->xcond_dims, "xcond_dims", &xcond_dims);

    char *c_ptr = (char *) raw_memory;

    ocp_qp_full_condensing_solver_memory *mem = (ocp_qp_full_condensing_solver_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_full_condensing_solver_memory);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    mem->cond_memory = ocp_qp_full_condensing_memory_assign(dims->xcond_dims, opts->xcond_opts, c_ptr);
    c_ptr += ocp_qp_full_condensing_memory_calculate_size(dims->xcond_dims, opts->xcond_opts);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    mem->solver_memory = qp_solver->memory_assign(qp_solver, xcond_dims, opts->qp_solver_opts, c_ptr);
    c_ptr += qp_solver->memory_calculate_size(qp_solver, xcond_dims, opts->qp_solver_opts);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    mem->qpd_in = dense_qp_in_assign(xcond_dims, c_ptr);
    c_ptr += dense_qp_in_calculate_size(xcond_dims);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    mem->qpd_out = dense_qp_out_assign(xcond_dims, c_ptr);
    c_ptr += dense_qp_out_calculate_size(xcond_dims);

    assert((char *) raw_memory + ocp_qp_full_condensing_solver_memory_calculate_size(config_, dims, opts_) == c_ptr);

    return mem;
}
#endif



/************************************************
 * workspace
 ************************************************/

#if 0
int ocp_qp_full_condensing_solver_workspace_calculate_size(void *config_, ocp_qp_xcond_solver_dims *dims,
                                                           void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
	ocp_qp_xcond_config *xcond = config->xcond;

//    ocp_qp_full_condensing_solver_opts *opts = (ocp_qp_full_condensing_solver_opts *) opts_;
    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;

    void *xcond_dims;
	xcond->dims_get(xcond, dims->xcond_dims, "xcond_dims", &xcond_dims);

    int size = sizeof(ocp_qp_full_condensing_solver_workspace);
    size += ocp_qp_full_condensing_workspace_calculate_size(dims, opts->xcond_opts);
    size += qp_solver->workspace_calculate_size(qp_solver, xcond_dims, opts->qp_solver_opts);

    return size;
}
#endif



static void cast_workspace(void *config_, ocp_qp_xcond_solver_dims *dims,
                           void *opts_,
                           ocp_qp_full_condensing_solver_memory *mem,
                           ocp_qp_full_condensing_solver_workspace *work)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
	ocp_qp_xcond_config *xcond = config->xcond;

    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;

    void *xcond_dims;
	xcond->dims_get(xcond, dims->xcond_dims, "xcond_dims", &xcond_dims);

    char *c_ptr = (char *) work;

    c_ptr += sizeof(ocp_qp_full_condensing_solver_workspace);

    work->cond_work = c_ptr;
    c_ptr += ocp_qp_full_condensing_workspace_calculate_size(dims, opts->xcond_opts);

    work->solver_workspace = c_ptr;
    c_ptr += qp_solver->workspace_calculate_size(qp_solver, xcond_dims, opts->qp_solver_opts);

    assert((char *) work + config->workspace_calculate_size(config_, dims, opts) >= c_ptr);
}



/************************************************
 * functions
 ************************************************/

int ocp_qp_full_condensing_solver(void *config_, ocp_qp_xcond_solver_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts_,
                                  void *mem_, void *work_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    qp_info *info = (qp_info *) qp_out->misc;
    acados_timer tot_timer, cond_timer;
    acados_tic(&tot_timer);

    // cast data structures
//    ocp_qp_full_condensing_solver_opts *opts = (ocp_qp_full_condensing_solver_opts *) opts_;
    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;
    ocp_qp_full_condensing_solver_memory *memory = (ocp_qp_full_condensing_solver_memory *) mem_;
    ocp_qp_full_condensing_solver_workspace *work =
        (ocp_qp_full_condensing_solver_workspace *) work_;

    // cast workspace
    cast_workspace(config_, dims, opts, memory, work);

    // condense
    acados_tic(&cond_timer);
    ocp_qp_full_condensing(qp_in, memory->qpd_in, opts->xcond_opts, memory->cond_memory,
                           work->cond_work);
    info->condensing_time = acados_toc(&cond_timer);

    // solve qp
    int solver_status =
        qp_solver->evaluate(qp_solver, memory->qpd_in, memory->qpd_out, opts->qp_solver_opts,
                            memory->solver_memory, work->solver_workspace);

    // expand
    acados_tic(&cond_timer);
    ocp_qp_full_expansion(memory->qpd_out, qp_out, opts->xcond_opts, memory->cond_memory,
                          work->cond_work);
    info->condensing_time += acados_toc(&cond_timer);

	// TODO something different to work with any QP type .....
	qp_info *info_mem = (qp_info *) memory->qpd_out->misc;
    info->total_time = acados_toc(&tot_timer);
    info->solve_QP_time = info_mem->solve_QP_time;
    info->interface_time = info_mem->interface_time;
    info->num_iter = info_mem->num_iter;
    info->t_computed = info_mem->t_computed;

    return solver_status;
}



void ocp_qp_full_condensing_solver_eval_sens(void *config_, ocp_qp_xcond_solver_dims *dims, ocp_qp_in *param_qp_in, ocp_qp_out *sens_qp_out,
                                     void *opts_, void *mem_, void *work_)
{
	printf("\nocp_qp_full_condensing_solver_eval_sens: not implemented yet\n");
	exit(1);
	return;
}



void ocp_qp_full_condensing_solver_config_initialize_default(void *config_)
{
    ocp_qp_xcond_solver_config *config = config_;

//    config->dims_set = &ocp_qp_dims_set;
//    config->opts_calculate_size = &ocp_qp_full_condensing_solver_opts_calculate_size;
//    config->opts_assign = &ocp_qp_full_condensing_solver_opts_assign;
//    config->opts_set = &ocp_qp_full_condensing_solver_opts_set;
//    config->opts_initialize_default = &ocp_qp_full_condensing_solver_opts_initialize_default;
//    config->opts_update = &ocp_qp_full_condensing_solver_opts_update;
    config->dims_calculate_size = &ocp_qp_xcond_solver_dims_calculate_size;
    config->dims_assign = &ocp_qp_xcond_solver_dims_assign;
    config->dims_set = &ocp_qp_xcond_solver_dims_set;
    config->opts_calculate_size = &ocp_qp_xcond_solver_opts_calculate_size;
    config->opts_assign = &ocp_qp_xcond_solver_opts_assign;
    config->opts_initialize_default = &ocp_qp_xcond_solver_opts_initialize_default;
    config->opts_update = &ocp_qp_xcond_solver_opts_update;
    config->opts_set = &ocp_qp_xcond_solver_opts_set;
    config->memory_calculate_size = &ocp_qp_xcond_solver_memory_calculate_size;
    config->memory_assign = &ocp_qp_xcond_solver_memory_assign;
    config->workspace_calculate_size = &ocp_qp_xcond_solver_workspace_calculate_size;

//    config->memory_calculate_size = &ocp_qp_full_condensing_solver_memory_calculate_size;
//    config->memory_assign = &ocp_qp_full_condensing_solver_memory_assign;
//    config->workspace_calculate_size = &ocp_qp_full_condensing_solver_workspace_calculate_size;
    config->evaluate = &ocp_qp_full_condensing_solver;
    config->eval_sens = &ocp_qp_full_condensing_solver_eval_sens;

	// initialize xcond
	ocp_qp_full_condensing_config_initialize_default(config->xcond);

    return;
}
