/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren, Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor, Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan, Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// external
#include <assert.h>
#include <string.h>
#include <stdlib.h>

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing.h" // TODO remove !!!
#include "acados/ocp_qp/ocp_qp_partial_condensing_solver.h" // TODO rename !!!
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"



/************************************************
 * opts
 ************************************************/

int ocp_qp_xcond_solver_opts_calculate_size(void *config_, ocp_qp_dims *dims)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
	ocp_qp_xcond_config *xcond = config->xcond;

    int size = 0;
    size += sizeof(ocp_qp_xcond_solver_opts);

    size += xcond->opts_calculate_size(dims);
    size += qp_solver->opts_calculate_size(qp_solver, dims);

    return size;
}



void *ocp_qp_xcond_solver_opts_assign(void *config_, ocp_qp_dims *dims, void *raw_memory)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
	ocp_qp_xcond_config *xcond = config->xcond;

    char *c_ptr = (char *) raw_memory;

    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) c_ptr;
    c_ptr += sizeof(ocp_qp_xcond_solver_opts);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    opts->xcond_opts = xcond->opts_assign(dims, c_ptr);
    c_ptr += xcond->opts_calculate_size(dims);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    opts->qp_solver_opts = qp_solver->opts_assign(qp_solver, dims, c_ptr);
    c_ptr += qp_solver->opts_calculate_size(qp_solver, dims);

    assert((char *) raw_memory +
               ocp_qp_xcond_solver_opts_calculate_size(config_, dims) ==
           c_ptr);

    return (void *) opts;
}



void ocp_qp_xcond_solver_opts_initialize_default(void *config_, ocp_qp_dims *dims,
                                                              void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
	ocp_qp_xcond_config *xcond = config->xcond;

    // xcond solver opts
    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;
    // xcond opts
    xcond->opts_initialize_default(dims, opts->xcond_opts);
    // qp solver opts
    qp_solver->opts_initialize_default(qp_solver, dims, opts->qp_solver_opts);
}



void ocp_qp_xcond_solver_opts_update(void *config_, ocp_qp_dims *dims, void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
	ocp_qp_xcond_config *xcond = config->xcond;

    // xcond solver opts
    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;
    // xcond opts
    xcond->opts_update(dims, opts->xcond_opts);
    // qp solver opts
    qp_solver->opts_update(qp_solver, dims, opts->qp_solver_opts);
}



void ocp_qp_xcond_solver_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
	ocp_qp_xcond_config *xcond = config->xcond;

	int ii;

	char module[MAX_STR_LEN];
	char *ptr_module = NULL;
	int module_length = 0;

	// extract module name
	char *char_ = strchr(field, '_');
	if(char_!=NULL)
	{
		module_length = char_-field;
		for(ii=0; ii<module_length; ii++)
			module[ii] = field[ii];
		module[module_length] = '\0'; // add end of string
		ptr_module = module;
	}

	if(!strcmp(ptr_module, "cond")) // pass options to condensing module
	{
		xcond->opts_set(opts->xcond_opts, field+module_length+1, value);
	}
	else // pass options to QP module
	{
		qp_solver->opts_set(qp_solver, opts->qp_solver_opts, field, value);
	}

	return;

}



/************************************************
 * memory
 ************************************************/

int ocp_qp_partial_condensing_solver_memory_calculate_size(void *config_, ocp_qp_dims *dims,
                                                           void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
	ocp_qp_xcond_config *xcond = config->xcond;

    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;

    int size = 0;
    size += sizeof(ocp_qp_partial_condensing_solver_memory);

    // set up dimesions of partially condensed qp
    void *xcond_dims;
	xcond->opts_get(opts->xcond_opts, "xcond_dims", &xcond_dims);
//    if (opts->xcond_opts->N2 < dims->N)
//    {
//        pcond_dims = opts->xcond_opts->pcond_dims;
//    }
//    else
//    {
//        pcond_dims = dims;
//    }

//    if (opts->xcond_opts->N2 < dims->N)
//    {
        size += ocp_qp_partial_condensing_memory_calculate_size(dims, opts->xcond_opts);
//    }

    size += qp_solver->memory_calculate_size(qp_solver, xcond_dims, opts->qp_solver_opts);

//    if (opts->xcond_opts->N2 < dims->N)
//    {
        size += ocp_qp_in_calculate_size(qp_solver, xcond_dims);
        size += ocp_qp_out_calculate_size(qp_solver, xcond_dims);
//    }

    return size;
}



void *ocp_qp_partial_condensing_solver_memory_assign(void *config_, ocp_qp_dims *dims, void *opts_,
                                                     void *raw_memory)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
	ocp_qp_xcond_config *xcond = config->xcond;

    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;

//    assert(opts->xcond_opts->N2 > 0 && "N2 must be positive!");
//    assert(opts->xcond_opts->N2 <= dims->N && "N2 cannot be bigger than N!");

    char *c_ptr = (char *) raw_memory;

    // set up dimesions of partially condensed qp
    void *xcond_dims;
	xcond->opts_get(opts->xcond_opts, "xcond_dims", &xcond_dims);
//    if (opts->xcond_opts->N2 < dims->N)
//    {
//        pcond_dims = opts->xcond_opts->pcond_dims;
//    }
//    else
//    {
//        pcond_dims = dims;
//    }

    ocp_qp_partial_condensing_solver_memory *mem =
        (ocp_qp_partial_condensing_solver_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_partial_condensing_solver_memory);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

//    if (opts->xcond_opts->N2 < dims->N)
//    {
        mem->pcond_memory =
            (ocp_qp_partial_condensing_memory *) ocp_qp_partial_condensing_memory_assign(
                dims, opts->xcond_opts, c_ptr);
        c_ptr += ocp_qp_partial_condensing_memory_calculate_size(dims, opts->xcond_opts);
//    }
//    else
//    {
//        mem->pcond_memory = NULL;
//    }

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    mem->solver_memory =
        qp_solver->memory_assign(qp_solver, xcond_dims, opts->qp_solver_opts, c_ptr);
    c_ptr += qp_solver->memory_calculate_size(qp_solver, xcond_dims, opts->qp_solver_opts);

//    if (opts->xcond_opts->N2 < dims->N)
//    {
        mem->pcond_qp_in = ocp_qp_in_assign(qp_solver, xcond_dims, c_ptr);
        c_ptr += ocp_qp_in_calculate_size(qp_solver, xcond_dims);
//    }
//    else
//    {
//        mem->pcond_qp_in = NULL;
//    }

//    if (opts->xcond_opts->N2 < dims->N)
//    {
        mem->pcond_qp_out = ocp_qp_out_assign(qp_solver, xcond_dims, c_ptr);
        c_ptr += ocp_qp_out_calculate_size(qp_solver, xcond_dims);
//    }
//    else
//    {
//        mem->pcond_qp_out = NULL;
//    }

    assert((char *) raw_memory +
               ocp_qp_partial_condensing_solver_memory_calculate_size(config_, dims, opts_) ==
           c_ptr);

    return mem;
}



/************************************************
 * workspace
 ************************************************/

int ocp_qp_xcond_solver_workspace_calculate_size(void *config_, ocp_qp_dims *dims,
                                                              void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
	ocp_qp_xcond_config *xcond = config->xcond;

    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;

    int size = sizeof(ocp_qp_xcond_solver_workspace);

    size += xcond->workspace_calculate_size(dims, opts->xcond_opts);

    // set up dimesions of condensed qp
    void *xcond_dims;
	xcond->opts_get(opts->xcond_opts, "xcond_dims", &xcond_dims);
//    if (opts->xcond_opts->N2 < dims->N)
//    {
//        pcond_dims = opts->xcond_opts->pcond_dims;
//    }
//    else
//    {
//        pcond_dims = dims;
//    }

    size += qp_solver->workspace_calculate_size(qp_solver, xcond_dims, opts->qp_solver_opts);

    return size;
}



static void cast_workspace(void *config_, ocp_qp_dims *dims,
                           ocp_qp_xcond_solver_opts *opts,
                           ocp_qp_partial_condensing_solver_memory *mem,
                           ocp_qp_xcond_solver_workspace *work)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
	ocp_qp_xcond_config *xcond = config->xcond;

    // set up dimesions of  condensed qp
    void *xcond_dims;
	xcond->opts_get(opts->xcond_opts, "xcond_dims", &xcond_dims);
//    if (opts->xcond_opts->N2 < dims->N)
//    {
//        pcond_dims = opts->xcond_opts->pcond_dims;
//    }
//    else
//    {
//        pcond_dims = dims;
//    }

    char *c_ptr = (char *) work;

    c_ptr += sizeof(ocp_qp_xcond_solver_workspace);

    work->xcond_work = c_ptr;
    c_ptr += xcond->workspace_calculate_size(dims, opts->xcond_opts);

    work->qp_solver_work = c_ptr;
    c_ptr += qp_solver->workspace_calculate_size(qp_solver, xcond_dims, opts->qp_solver_opts);
}



/************************************************
 * functions
 ************************************************/

int ocp_qp_partial_condensing_solver(void *config_, ocp_qp_in *qp_in, ocp_qp_out *qp_out,
                                     void *opts_, void *mem_, void *work_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
	ocp_qp_xcond_config *xcond = config->xcond;

    ocp_qp_info *info = (ocp_qp_info *) qp_out->misc;
    acados_timer tot_timer, cond_timer;
    acados_tic(&tot_timer);

    // cast data structures
    ocp_qp_xcond_solver_opts *opts = opts_;
    ocp_qp_partial_condensing_solver_memory *memory = mem_;
    ocp_qp_xcond_solver_workspace *work = work_;

    // cast workspace
    cast_workspace(config_, qp_in->dim, opts, memory, work);

    int solver_status;
	int tmp_status;

//    if (opts->xcond_opts->N2 < qp_in->dim->N)
//    {  // condensing
        acados_tic(&cond_timer);
        tmp_status = xcond->condensing(qp_in, memory->pcond_qp_in, opts->xcond_opts, memory->pcond_memory, work->xcond_work);
        info->condensing_time = acados_toc(&cond_timer);
//    }
//    else
//    {
//        memory->pcond_qp_in = qp_in;
//        memory->pcond_qp_out = qp_out;
//        info->condensing_time = 0;
//    }

    // solve qp
	solver_status = qp_solver->evaluate(qp_solver, memory->pcond_qp_in, memory->pcond_qp_out,
        		opts->qp_solver_opts, memory->solver_memory, work->qp_solver_work);

    // expand
//    if (opts->xcond_opts->N2 < qp_in->dim->N)
//    {
        acados_tic(&cond_timer);

        tmp_status = xcond->expansion(memory->pcond_qp_out, qp_out, opts->xcond_opts,
				memory->pcond_memory, work->xcond_work);

        info->condensing_time += acados_toc(&cond_timer);
//    }

    info->total_time = acados_toc(&tot_timer);
    info->solve_QP_time = ((ocp_qp_info *) (memory->pcond_qp_out->misc))->solve_QP_time;
    info->interface_time = ((ocp_qp_info *) (memory->pcond_qp_out->misc))->interface_time;
    info->num_iter = ((ocp_qp_info *) (memory->pcond_qp_out->misc))->num_iter;
    info->t_computed = ((ocp_qp_info *) (memory->pcond_qp_out->misc))->t_computed;

    return solver_status;
}



void ocp_qp_partial_condensing_solver_eval_sens(void *config_, ocp_qp_in *param_qp_in, ocp_qp_out *sens_qp_out,
		void *opts_, void *mem_, void *work_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;

    // cast data structures
    ocp_qp_xcond_solver_opts *opts = opts_;
    ocp_qp_partial_condensing_solver_memory *memory = mem_;
    ocp_qp_xcond_solver_workspace *work = work_;

    // cast workspace
    cast_workspace(config_, param_qp_in->dim, opts, memory, work);

	int tmp_status;

	// condensing
//    if (opts->xcond_opts->N2 < param_qp_in->dim->N)
//    {
        tmp_status = ocp_qp_partial_condensing(param_qp_in, memory->pcond_qp_in, opts->xcond_opts,
        		memory->pcond_memory, work->xcond_work);
//    }
//    else
//    {
//        memory->pcond_qp_in = param_qp_in;
//        memory->pcond_qp_out = sens_qp_out;
//    }

	// eval sensitivity
//	solver_status = qp_solver->evaluate(qp_solver, memory->pcond_qp_in, memory->pcond_qp_out,
//			opts->qp_solver_opts, memory->solver_memory, work->qp_solver_work);

	printf("\nocp_qp_partial_condensing_solver_eval_sens: not implemented yet\n");
	exit(1);

	return;

}



void ocp_qp_partial_condensing_solver_config_initialize_default(void *config_)
{
    ocp_qp_xcond_solver_config *config = config_;

    config->dims_set = &ocp_qp_dims_set;
    config->opts_calculate_size = &ocp_qp_xcond_solver_opts_calculate_size;
    config->opts_assign = &ocp_qp_xcond_solver_opts_assign;
    config->opts_initialize_default = &ocp_qp_xcond_solver_opts_initialize_default;
    config->opts_update = &ocp_qp_xcond_solver_opts_update;
    config->opts_set = &ocp_qp_xcond_solver_opts_set;
    config->memory_calculate_size = &ocp_qp_partial_condensing_solver_memory_calculate_size;
    config->memory_assign = &ocp_qp_partial_condensing_solver_memory_assign;
    config->workspace_calculate_size = &ocp_qp_xcond_solver_workspace_calculate_size;
    config->evaluate = &ocp_qp_partial_condensing_solver;
    config->eval_sens = &ocp_qp_partial_condensing_solver_eval_sens;

    return;
}
