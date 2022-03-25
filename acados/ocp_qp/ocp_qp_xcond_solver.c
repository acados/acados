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
#include <string.h>
#include <stdlib.h>

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_xcond_solver.h"
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"



/************************************************
 * config
 ************************************************/

acados_size_t ocp_qp_xcond_solver_config_calculate_size()
{
    acados_size_t size = 0;

    size += sizeof(ocp_qp_xcond_solver_config);

    size += ocp_qp_solver_config_calculate_size();  // qp_solver
    size += ocp_qp_condensing_config_calculate_size();  // xcond

    return size;
}



ocp_qp_xcond_solver_config *ocp_qp_xcond_solver_config_assign(void *raw_memory)
{
    char *c_ptr = raw_memory;

    ocp_qp_xcond_solver_config *config = (ocp_qp_xcond_solver_config *) c_ptr;
    c_ptr += sizeof(ocp_qp_xcond_solver_config);

    config->qp_solver = ocp_qp_solver_config_assign(c_ptr);
    c_ptr += ocp_qp_solver_config_calculate_size();

    config->xcond = ocp_qp_condensing_config_assign(c_ptr);
    c_ptr += ocp_qp_condensing_config_calculate_size();

    return config;
}



/************************************************
 * dims
 ************************************************/

acados_size_t ocp_qp_xcond_solver_dims_calculate_size(void *config_, int N)
{
    ocp_qp_xcond_solver_config *config = config_;

    acados_size_t size = sizeof(ocp_qp_xcond_solver_dims);

    // orig_dims
    size += ocp_qp_dims_calculate_size(N);

    // xcond_dims
    size += config->xcond->dims_calculate_size(config->xcond, N);

    return size;
}



ocp_qp_xcond_solver_dims *ocp_qp_xcond_solver_dims_assign(void *config_, int N, void *raw_memory)
{
    ocp_qp_xcond_solver_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    ocp_qp_xcond_solver_dims *dims = (ocp_qp_xcond_solver_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_xcond_solver_dims);

    // orig_dims
    dims->orig_dims = ocp_qp_dims_assign(N, c_ptr);
    c_ptr += ocp_qp_dims_calculate_size(N);

    // xcond_dims
    dims->xcond_dims = config->xcond->dims_assign(config->xcond, N, c_ptr);
    c_ptr += config->xcond->dims_calculate_size(config->xcond, N);

    assert((char *) raw_memory + ocp_qp_xcond_solver_dims_calculate_size(config_, N) == c_ptr);

    return dims;
}



void ocp_qp_xcond_solver_dims_set_(void *config_, ocp_qp_xcond_solver_dims *dims,
                                  int stage, const char *field, int* value)
{
    ocp_qp_xcond_solver_config *config = config_;

    // orig_dims
    ocp_qp_dims_set(config_, dims->orig_dims, stage, field, value);

    // xcond_dims
    config->xcond->dims_set(config->xcond, dims->xcond_dims, stage, field, value);

    return;
}



/************************************************
 * opts
 ************************************************/

acados_size_t ocp_qp_xcond_solver_opts_calculate_size(void *config_, ocp_qp_xcond_solver_dims *dims)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
    ocp_qp_xcond_config *xcond = config->xcond;

    void *xcond_qp_dims;
    xcond->dims_get(xcond, dims->xcond_dims, "xcond_dims", &xcond_qp_dims);

    acados_size_t size = 0;

    size += sizeof(ocp_qp_xcond_solver_opts);

    size += xcond->opts_calculate_size(dims->xcond_dims);

    size += qp_solver->opts_calculate_size(qp_solver, xcond_qp_dims);

    return size;
}



void *ocp_qp_xcond_solver_opts_assign(void *config_, ocp_qp_xcond_solver_dims *dims, void *raw_memory)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
    ocp_qp_xcond_config *xcond = config->xcond;

    void *xcond_qp_dims;
    xcond->dims_get(xcond, dims->xcond_dims, "xcond_dims", &xcond_qp_dims);

    char *c_ptr = (char *) raw_memory;

    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) c_ptr;
    c_ptr += sizeof(ocp_qp_xcond_solver_opts);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    opts->xcond_opts = xcond->opts_assign(dims->xcond_dims, c_ptr);
    c_ptr += xcond->opts_calculate_size(dims->xcond_dims);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    opts->qp_solver_opts = qp_solver->opts_assign(qp_solver, xcond_qp_dims, c_ptr);
    c_ptr += qp_solver->opts_calculate_size(qp_solver, xcond_qp_dims);

    assert((char *) raw_memory + ocp_qp_xcond_solver_opts_calculate_size(config_, dims) == c_ptr);

    return (void *) opts;
}



void ocp_qp_xcond_solver_opts_initialize_default(void *config_, ocp_qp_xcond_solver_dims *dims, void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
    ocp_qp_xcond_config *xcond = config->xcond;

    void *xcond_qp_dims;
    xcond->dims_get(xcond, dims->xcond_dims, "xcond_dims", &xcond_qp_dims);

    // xcond solver opts
    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;
    // xcond opts
    xcond->opts_initialize_default(dims->xcond_dims, opts->xcond_opts);
    // qp solver opts
    qp_solver->opts_initialize_default(qp_solver, xcond_qp_dims, opts->qp_solver_opts);
}



void ocp_qp_xcond_solver_opts_update(void *config_, ocp_qp_xcond_solver_dims *dims, void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
    ocp_qp_xcond_config *xcond = config->xcond;

    void *xcond_qp_dims;
    xcond->dims_get(xcond, dims->xcond_dims, "xcond_dims", &xcond_qp_dims);

    // xcond solver opts
    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;
    // xcond opts
    xcond->opts_update(dims->xcond_dims, opts->xcond_opts);
    // qp solver opts
    qp_solver->opts_update(qp_solver, xcond_qp_dims, opts->qp_solver_opts);
}



void ocp_qp_xcond_solver_opts_set_(void *config_, void *opts_, const char *field, void* value)
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

    if( ptr_module!=NULL && (!strcmp(ptr_module, "cond")) ) // pass options to condensing module // TODO rename xcond ???
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

acados_size_t ocp_qp_xcond_solver_memory_calculate_size(void *config_, ocp_qp_xcond_solver_dims *dims, void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
    ocp_qp_xcond_config *xcond = config->xcond;

    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;

    acados_size_t size = 0;
    size += sizeof(ocp_qp_xcond_solver_memory);

    // set up dimesions of partially condensed qp
    void *xcond_qp_dims;
    xcond->dims_get(xcond, dims->xcond_dims, "xcond_dims", &xcond_qp_dims);

    size += xcond->memory_calculate_size(dims->xcond_dims, opts->xcond_opts);

    size += qp_solver->memory_calculate_size(qp_solver, xcond_qp_dims, opts->qp_solver_opts);

    return size;
}



void *ocp_qp_xcond_solver_memory_assign(void *config_, ocp_qp_xcond_solver_dims *dims, void *opts_,
                                                     void *raw_memory)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
    ocp_qp_xcond_config *xcond = config->xcond;

    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;

    char *c_ptr = (char *) raw_memory;

    // set up dimesions of partially condensed qp
    void *xcond_qp_dims;
    xcond->dims_get(xcond, dims->xcond_dims, "xcond_dims", &xcond_qp_dims);

    ocp_qp_xcond_solver_memory *mem = (ocp_qp_xcond_solver_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_xcond_solver_memory);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    mem->xcond_memory = xcond->memory_assign(dims->xcond_dims, opts->xcond_opts, c_ptr);
    c_ptr += xcond->memory_calculate_size(dims->xcond_dims, opts->xcond_opts);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    mem->solver_memory = qp_solver->memory_assign(qp_solver, xcond_qp_dims, opts->qp_solver_opts, c_ptr);
    c_ptr += qp_solver->memory_calculate_size(qp_solver, xcond_qp_dims, opts->qp_solver_opts);

    xcond->memory_get(xcond, mem->xcond_memory, "xcond_qp_in", &mem->xcond_qp_in);
    xcond->memory_get(xcond, mem->xcond_memory, "xcond_qp_out", &mem->xcond_qp_out);

    assert((char *) raw_memory + ocp_qp_xcond_solver_memory_calculate_size(config_, dims, opts_) >= c_ptr);

    return mem;
}



void ocp_qp_xcond_solver_memory_reset(void *config_, ocp_qp_xcond_solver_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out,
                                     void *opts_, void *mem_, void *work_)
{
    // cast data structures
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
    ocp_qp_xcond_config *xcond = config->xcond;
    ocp_qp_xcond_solver_opts *opts = opts_;
    ocp_qp_xcond_solver_memory *mem = mem_;
    // ocp_qp_xcond_solver_workspace *work = work_;

    void *xcond_qp_dims;
    xcond->dims_get(xcond, dims->xcond_dims, "xcond_dims", &xcond_qp_dims);

    mem->solver_memory = qp_solver->memory_assign(qp_solver, xcond_qp_dims, opts->qp_solver_opts, mem->solver_memory);

    return;
}



void ocp_qp_xcond_solver_memory_get(void *config_, void *mem_, const char *field, void* value)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
    ocp_qp_xcond_config *xcond = config->xcond;

    ocp_qp_xcond_solver_memory *mem = mem_;

    // TODO extract module name as for opts_set

    if (!strcmp(field, "time_qp_solver_call"))
    {
        qp_solver->memory_get(qp_solver, mem->solver_memory, field, value);
    }
    else if (!strcmp(field, "iter"))
    {
        qp_solver->memory_get(qp_solver, mem->solver_memory, field, value);
    }
    else if (!strcmp(field, "status"))
    {
        qp_solver->memory_get(qp_solver, mem->solver_memory, field, value);
    }
    else if (!strcmp(field, "time_qp_xcond"))
    {
        xcond->memory_get(xcond, mem->xcond_memory, field, value);
    }
    else
    {
        printf("\nerror: ocp_qp_xcond_solver_memory_get: field %s not available\n", field);
        exit(1);
    }

    return;

}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_qp_xcond_solver_workspace_calculate_size(void *config_, ocp_qp_xcond_solver_dims *dims, void *opts_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
    ocp_qp_xcond_config *xcond = config->xcond;

    ocp_qp_xcond_solver_opts *opts = (ocp_qp_xcond_solver_opts *) opts_;

    acados_size_t size = sizeof(ocp_qp_xcond_solver_workspace);

//    size += xcond->workspace_calculate_size(dims->orig_dims, opts->xcond_opts);
    size += xcond->workspace_calculate_size(dims->xcond_dims, opts->xcond_opts);

    // set up dimesions of condensed qp
    void *xcond_qp_dims;
    xcond->dims_get(xcond, dims->xcond_dims, "xcond_dims", &xcond_qp_dims);

    size += qp_solver->workspace_calculate_size(qp_solver, xcond_qp_dims, opts->qp_solver_opts);

    return size;
}



static void cast_workspace(void *config_, ocp_qp_xcond_solver_dims *dims,
                           ocp_qp_xcond_solver_opts *opts,
                           ocp_qp_xcond_solver_memory *mem,
                           ocp_qp_xcond_solver_workspace *work)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
    ocp_qp_xcond_config *xcond = config->xcond;

    // set up dimesions of  condensed qp
    void *xcond_qp_dims;
    xcond->dims_get(xcond, dims->xcond_dims, "xcond_dims", &xcond_qp_dims);

    char *c_ptr = (char *) work;

    c_ptr += sizeof(ocp_qp_xcond_solver_workspace);

    work->xcond_work = c_ptr;
    c_ptr += xcond->workspace_calculate_size(dims->xcond_dims, opts->xcond_opts);

    work->qp_solver_work = c_ptr;
    c_ptr += qp_solver->workspace_calculate_size(qp_solver, xcond_qp_dims, opts->qp_solver_opts);

    assert((char *) work + config->workspace_calculate_size(config_, dims, opts) >= c_ptr);
}



/************************************************
 * functions
 ************************************************/

int ocp_qp_xcond_solver(void *config_, ocp_qp_xcond_solver_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out,
                                     void *opts_, void *mem_, void *work_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
    ocp_qp_xcond_config *xcond = config->xcond;

    qp_info *info = (qp_info *) qp_out->misc;
    acados_timer tot_timer, cond_timer;
    acados_tic(&tot_timer);

    // cast data structures
    ocp_qp_xcond_solver_opts *opts = opts_;
    ocp_qp_xcond_solver_memory *memory = mem_;
    ocp_qp_xcond_solver_workspace *work = work_;

    // cast workspace
    cast_workspace(config_, dims, opts, memory, work);

    int solver_status = ACADOS_SUCCESS;

    // condensing
    acados_tic(&cond_timer);
    xcond->condensing(qp_in, memory->xcond_qp_in, opts->xcond_opts, memory->xcond_memory, work->xcond_work);
    info->condensing_time = acados_toc(&cond_timer);

    // solve qp
    solver_status = qp_solver->evaluate(qp_solver, memory->xcond_qp_in, memory->xcond_qp_out,
                                opts->qp_solver_opts, memory->solver_memory, work->qp_solver_work);

    // expansion
    acados_tic(&cond_timer);
    xcond->expansion(memory->xcond_qp_out, qp_out, opts->xcond_opts, memory->xcond_memory, work->xcond_work);
    info->condensing_time += acados_toc(&cond_timer);

    // output qp info
    qp_info *info_mem;
    xcond->memory_get(xcond, memory->xcond_memory, "qp_out_info", &info_mem);

    info->total_time = acados_toc(&tot_timer);
    info->solve_QP_time = info_mem->solve_QP_time;
    info->interface_time = info_mem->interface_time;
    info->num_iter = info_mem->num_iter;
    info->t_computed = info_mem->t_computed;

    return solver_status;
}



void ocp_qp_xcond_solver_eval_sens(void *config_, ocp_qp_xcond_solver_dims *dims, ocp_qp_in *param_qp_in, ocp_qp_out *sens_qp_out,
        void *opts_, void *mem_, void *work_)
{
    ocp_qp_xcond_solver_config *config = config_;
    qp_solver_config *qp_solver = config->qp_solver;
    ocp_qp_xcond_config *xcond = config->xcond;

//    qp_info *info = (qp_info *) qp_out->misc;
//    acados_timer tot_timer, cond_timer;
//    acados_tic(&tot_timer);

    // cast data structures
    ocp_qp_xcond_solver_opts *opts = opts_;
    ocp_qp_xcond_solver_memory *memory = mem_;
    ocp_qp_xcond_solver_workspace *work = work_;

    // cast workspace
    cast_workspace(config_, dims, opts, memory, work);


    // condensing
//    acados_tic(&cond_timer);
    xcond->condensing_rhs(param_qp_in, memory->xcond_qp_in, opts->xcond_opts, memory->xcond_memory, work->xcond_work);
//    info->condensing_time = acados_toc(&cond_timer);

    // qp evaluate sensitivity
    qp_solver->eval_sens(qp_solver, memory->xcond_qp_in, memory->xcond_qp_out, opts->qp_solver_opts, memory->solver_memory, work->qp_solver_work);

    // expansion
//    acados_tic(&cond_timer);
    xcond->expansion(memory->xcond_qp_out, sens_qp_out, opts->xcond_opts, memory->xcond_memory, work->xcond_work);
//    info->condensing_time += acados_toc(&cond_timer);

    // output qp info
//    qp_info *info_mem;
//    xcond->memory_get(xcond, memory->xcond_memory, "qp_out_info", &info_mem);

//    info->total_time = acados_toc(&tot_timer);
//    info->solve_QP_time = info_mem->solve_QP_time;
//    info->interface_time = info_mem->interface_time;
//    info->num_iter = info_mem->num_iter;
//    info->t_computed = info_mem->t_computed;

    return;

}



void ocp_qp_xcond_solver_config_initialize_default(void *config_)
{
    ocp_qp_xcond_solver_config *config = config_;

    config->dims_calculate_size = &ocp_qp_xcond_solver_dims_calculate_size;
    config->dims_assign = &ocp_qp_xcond_solver_dims_assign;
    config->dims_set = &ocp_qp_xcond_solver_dims_set_;
    config->opts_calculate_size = &ocp_qp_xcond_solver_opts_calculate_size;
    config->opts_assign = &ocp_qp_xcond_solver_opts_assign;
    config->opts_initialize_default = &ocp_qp_xcond_solver_opts_initialize_default;
    config->opts_update = &ocp_qp_xcond_solver_opts_update;
    config->opts_set = &ocp_qp_xcond_solver_opts_set_;
    config->memory_calculate_size = &ocp_qp_xcond_solver_memory_calculate_size;
    config->memory_assign = &ocp_qp_xcond_solver_memory_assign;
    config->memory_get = &ocp_qp_xcond_solver_memory_get;
    config->memory_reset = &ocp_qp_xcond_solver_memory_reset;
    config->workspace_calculate_size = &ocp_qp_xcond_solver_workspace_calculate_size;
    config->evaluate = &ocp_qp_xcond_solver;
    config->eval_sens = &ocp_qp_xcond_solver_eval_sens;

    return;
}

