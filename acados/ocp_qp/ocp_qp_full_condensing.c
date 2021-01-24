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
#include <stdlib.h>
#include <assert.h>
#include <string.h>
// hpipm
#include "hpipm/include/hpipm_d_cond.h"
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_ipm.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp_red.h"
// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_full_condensing.h"
#include "acados/utils/mem.h"
#include "acados/utils/types.h"
#include "acados/utils/timing.h"



/************************************************
 * dims
 ************************************************/

acados_size_t ocp_qp_full_condensing_dims_calculate_size(void *config, int N)
{
    acados_size_t size = 0;

    size += sizeof(ocp_qp_full_condensing_dims);

    // orig_dims, red_dims
    size += 2*sizeof(ocp_qp_dims);
    size += 2*d_ocp_qp_dim_memsize(N);

    // fcond_dims
    size += sizeof(dense_qp_dims);

    size += 1*8;

    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_qp_full_condensing_dims_assign(void *config, int N, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // dims
    ocp_qp_full_condensing_dims *dims = (ocp_qp_full_condensing_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_full_condensing_dims);

    // orig_dims
    dims->orig_dims = (ocp_qp_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_dims);
    // red_dims
    dims->red_dims = (ocp_qp_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_dims);
    // fcond_dims
    dims->fcond_dims = (dense_qp_dims *) c_ptr;
    c_ptr += sizeof(dense_qp_dims);

    align_char_to(8, &c_ptr);

    // orig_dims
    d_ocp_qp_dim_create(N, dims->orig_dims, c_ptr);
    c_ptr += d_ocp_qp_dim_memsize(N);
    // red_dims
    d_ocp_qp_dim_create(N, dims->red_dims, c_ptr);
    c_ptr += d_ocp_qp_dim_memsize(N);

    assert((char *) raw_memory + ocp_qp_full_condensing_dims_calculate_size(config, N) >= c_ptr);

    return dims;
}



void ocp_qp_full_condensing_dims_set(void *config_, void *dims_, int stage, const char *field, int *value)
{
    ocp_qp_full_condensing_dims *dims = dims_;

    ocp_qp_dims_set(config_, dims->orig_dims, stage, field, value);

    // TODO later in mem do the fcond_dims

    return;
}



void ocp_qp_full_condensing_dims_get(void *config_, void *dims_, const char *field, void* value)
{

    ocp_qp_full_condensing_dims *dims = dims_;

    if(!strcmp(field, "xcond_dims"))
    {
        dense_qp_dims **ptr = value;
        *ptr = dims->fcond_dims;
    }
    else
    {
        printf("\nerror: ocp_qp_full_condensing_dims_get: field %s not available\n", field);
        exit(1);
    }

    return;

}



/************************************************
 * opts
 ************************************************/

acados_size_t ocp_qp_full_condensing_opts_calculate_size(void *dims_)
{
    ocp_qp_full_condensing_dims *dims = dims_;

    // populate dimensions of reduced qp
    d_ocp_qp_dim_reduce_eq_dof(dims->orig_dims, dims->red_dims);
    // populate dimensions of new dense_qp
//    d_cond_qp_compute_dim(dims->orig_dims, dims->fcond_dims);
    d_cond_qp_compute_dim(dims->red_dims, dims->fcond_dims);
//d_ocp_qp_dim_print(dims->orig_dims);
//d_ocp_qp_dim_print(dims->red_dims);
//exit(1);

    acados_size_t size = 0;

    size += sizeof(ocp_qp_full_condensing_opts);

    // hpipm_cond_opts
    size += sizeof(struct d_cond_qp_arg);
    size += d_cond_qp_arg_memsize();
    // hpipm_red_opts
    size += sizeof(struct d_ocp_qp_reduce_eq_dof_arg);
    size += d_ocp_qp_reduce_eq_dof_arg_memsize();

    // fcond_dims
//    size += sizeof(dense_qp_dims);

    //
    size += 1*8;
    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_qp_full_condensing_opts_assign(void *dims_, void *raw_memory)
{
    // ocp_qp_full_condensing_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    // opts
    ocp_qp_full_condensing_opts *opts = (ocp_qp_full_condensing_opts *) c_ptr;
    c_ptr += sizeof(ocp_qp_full_condensing_opts);

    // fcond_dims
//    opts->fcond_dims = (dense_qp_dims *) c_ptr;
//    c_ptr += sizeof(dense_qp_dims);

    // hpipm_cond_opts
    opts->hpipm_cond_opts = (struct d_cond_qp_arg *) c_ptr;
    c_ptr += sizeof(struct d_cond_qp_arg);
    // hpipm_red_opts
    opts->hpipm_red_opts = (struct d_ocp_qp_reduce_eq_dof_arg *) c_ptr;
    c_ptr += sizeof(struct d_ocp_qp_reduce_eq_dof_arg);

    align_char_to(8, &c_ptr);

    // hpipm_cond_opts
    d_cond_qp_arg_create(opts->hpipm_cond_opts, c_ptr);
    c_ptr += opts->hpipm_cond_opts->memsize;
    // hpipm_red_opts
    d_ocp_qp_reduce_eq_dof_arg_create(opts->hpipm_red_opts, c_ptr);
    c_ptr += opts->hpipm_red_opts->memsize;

    assert((char *) raw_memory + ocp_qp_full_condensing_opts_calculate_size(dims_) >= c_ptr);

    return opts;
}



void ocp_qp_full_condensing_opts_initialize_default(void *dims_, void *opts_)
{
    ocp_qp_full_condensing_opts *opts = opts_;

    // condense both Hessian and gradient by default
    opts->cond_hess = 1;
    // expand only primal solution (linear MPC, Gauss-Newton)
    opts->expand_dual_sol = 1;

    // hpipm_cond_opts
    d_cond_qp_arg_set_default(opts->hpipm_cond_opts);

    // hpipm_red_opts
    d_ocp_qp_reduce_eq_dof_arg_set_default(opts->hpipm_red_opts);
    d_ocp_qp_reduce_eq_dof_arg_set_alias_unchanged(opts->hpipm_red_opts, 1);
    d_ocp_qp_reduce_eq_dof_arg_set_comp_dual_sol_eq(opts->hpipm_red_opts, 1);
    d_ocp_qp_reduce_eq_dof_arg_set_comp_dual_sol_ineq(opts->hpipm_red_opts, 1);

    opts->mem_qp_in = 1;

    return;
}



void ocp_qp_full_condensing_opts_update(void *dims_, void *opts_)
{
    ocp_qp_full_condensing_opts *opts = opts_;

    // hpipm_cond_opts
    d_cond_qp_arg_set_ric_alg(opts->ric_alg, opts->hpipm_cond_opts);

    return;
}



void ocp_qp_full_condensing_opts_set(void *opts_, const char *field, void* value)
{

    ocp_qp_full_condensing_opts *opts = opts_;

    if(!strcmp(field, "ric_alg"))
    {
        int *tmp_ptr = value;
        opts->ric_alg = *tmp_ptr;
    }
    else if(!strcmp(field, "hess"))
    {
        int *tmp_ptr = value;
        opts->cond_hess = *tmp_ptr;
    }
    else if(!strcmp(field, "dual_sol"))
    {
        int *tmp_ptr = value;
        opts->expand_dual_sol = *tmp_ptr;
        d_ocp_qp_reduce_eq_dof_arg_set_comp_dual_sol_eq(opts->hpipm_red_opts, *tmp_ptr);
        d_ocp_qp_reduce_eq_dof_arg_set_comp_dual_sol_ineq(opts->hpipm_red_opts, *tmp_ptr);
    }
    else
    {
        printf("\nerror: field %s not available in ocp_qp_full_condensing_opts_set\n", field);
        exit(1);
    }

    return;

}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_qp_full_condensing_memory_calculate_size(void *dims_, void *opts_)
{
    ocp_qp_full_condensing_dims *dims = dims_;
    ocp_qp_full_condensing_opts *opts = opts_;

    // TODO needed ???
    // populate dimensions of reduced qp
//    d_ocp_qp_dim_reduce_eq_dof(dims->orig_dims, dims->red_dims);
    // populate dimensions of new dense_qp
//    d_cond_qp_compute_dim(dims->orig_dims, dims->fcond_dims);
//    d_cond_qp_compute_dim(dims->red_dims, dims->fcond_dims);

    acados_size_t size = 0;

    size += sizeof(ocp_qp_full_condensing_memory);

    size += dense_qp_in_calculate_size(dims->fcond_dims);

    size += dense_qp_out_calculate_size(dims->fcond_dims);

    size += ocp_qp_in_calculate_size(dims->red_dims);

    size += ocp_qp_out_calculate_size(dims->red_dims);

    size += sizeof(struct d_cond_qp_ws);
    size += d_cond_qp_ws_memsize(dims->red_dims, opts->hpipm_cond_opts);

    size += sizeof(struct d_ocp_qp_reduce_eq_dof_ws);
    size += d_ocp_qp_reduce_eq_dof_ws_memsize(dims->orig_dims);

    size += 2*8;

    return size;
}



void *ocp_qp_full_condensing_memory_assign(void *dims_, void *opts_, void *raw_memory)
{
    ocp_qp_full_condensing_dims *dims = dims_;
    ocp_qp_full_condensing_opts *opts = opts_;

    char *c_ptr = (char *) raw_memory;

    // initial alignment
    align_char_to(8, &c_ptr);

    ocp_qp_full_condensing_memory *mem = (ocp_qp_full_condensing_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_full_condensing_memory);

    align_char_to(8, &c_ptr);

    // hpipm_cond_work struct
    mem->hpipm_cond_work = (struct d_cond_qp_ws *) c_ptr;
    c_ptr += sizeof(struct d_cond_qp_ws);
    // hpipm_red_work struct
    mem->hpipm_red_work = (struct d_ocp_qp_reduce_eq_dof_ws *) c_ptr;
    c_ptr += sizeof(struct d_ocp_qp_reduce_eq_dof_ws);

//    align_char_to(8, &c_ptr);
//    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    // hpipm_cond_work
    d_cond_qp_ws_create(dims->red_dims, opts->hpipm_cond_opts, mem->hpipm_cond_work, c_ptr);
    c_ptr += mem->hpipm_cond_work->memsize;
    // hpipm_red_work
    d_ocp_qp_reduce_eq_dof_ws_create(dims->orig_dims, mem->hpipm_red_work, c_ptr);
    c_ptr += mem->hpipm_red_work->memsize;

    mem->fcond_qp_in = dense_qp_in_assign(dims->fcond_dims, c_ptr);
    c_ptr += dense_qp_in_calculate_size(dims->fcond_dims);

    mem->fcond_qp_out = dense_qp_out_assign(dims->fcond_dims, c_ptr);
    c_ptr += dense_qp_out_calculate_size(dims->fcond_dims);

    mem->red_qp = ocp_qp_in_assign(dims->red_dims, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(dims->red_dims);

    mem->red_sol = ocp_qp_out_assign(dims->red_dims, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(dims->red_dims);

    mem->qp_out_info = (qp_info *) mem->fcond_qp_out->misc;

    assert((char *) raw_memory + ocp_qp_full_condensing_memory_calculate_size(dims, opts) >= c_ptr);

    return mem;
}



void ocp_qp_full_condensing_memory_get(void *config_, void *mem_, const char *field, void* value)
{
    ocp_qp_full_condensing_memory *mem = mem_;

    if(!strcmp(field, "xcond_qp_in"))
    {
        dense_qp_in **ptr = value;
        *ptr = mem->fcond_qp_in;
    }
    else if(!strcmp(field, "xcond_qp_out"))
    {
        dense_qp_out **ptr = value;
        *ptr = mem->fcond_qp_out;
    }
    else if(!strcmp(field, "qp_out_info"))
    {
        qp_info **ptr = value;
        *ptr = mem->qp_out_info;
    }
    else if (!strcmp(field, "time_qp_xcond"))
    {
        double *ptr = value;
        *ptr = mem->time_qp_xcond;
    }
    else
    {
        printf("\nerror: ocp_qp_full_condensing_memory_get: field %s not available\n", field);
        exit(1);
    }

    return;

}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_qp_full_condensing_workspace_calculate_size(void *dims_, void *opts_)
{

    return 0;

}



/************************************************
 * functions
 ************************************************/

int ocp_qp_full_condensing(void *qp_in_, void *fcond_qp_in_, void *opts_, void *mem_, void *work_)
{
    ocp_qp_in *qp_in = qp_in_;
    dense_qp_in *fcond_qp_in = fcond_qp_in_;
    ocp_qp_full_condensing_opts *opts = opts_;
    ocp_qp_full_condensing_memory *mem = mem_;

    acados_timer timer;

    // save pointer to ocp_qp_in in memory (needed for expansion)
    mem->ptr_qp_in = qp_in;

    // start timer
    acados_tic(&timer);

//d_ocp_qp_dim_print(qp_in->dim);
//d_ocp_qp_dim_print(mem->red_qp->dim);
    // reduce eq constr DOF
    d_ocp_qp_reduce_eq_dof(qp_in, mem->red_qp, opts->hpipm_red_opts, mem->hpipm_red_work);
//d_ocp_qp_print(qp_in->dim, qp_in);
//d_ocp_qp_print(mem->red_qp->dim, mem->red_qp);
//exit(1);

    // convert to dense qp structure
    if (opts->cond_hess == 0)
    {
        // condense gradient only
        d_cond_qp_cond_rhs(mem->red_qp, fcond_qp_in, opts->hpipm_cond_opts, mem->hpipm_cond_work);
    }
    else
    {
        // condense gradient and Hessian
        d_cond_qp_cond(mem->red_qp, fcond_qp_in, opts->hpipm_cond_opts, mem->hpipm_cond_work);
    }

    // stop timer
    mem->time_qp_xcond = acados_toc(&timer);

    return ACADOS_SUCCESS;
}



int ocp_qp_full_condensing_rhs(void *qp_in_, void *fcond_qp_in_, void *opts_, void *mem_, void *work_)
{
    ocp_qp_in *qp_in = qp_in_;
    dense_qp_in *fcond_qp_in = fcond_qp_in_;
    ocp_qp_full_condensing_opts *opts = opts_;
    ocp_qp_full_condensing_memory *mem = mem_;

    acados_timer timer;

    // start timer
    acados_tic(&timer);

    // save pointer to ocp_qp_in in memory (needed for expansion)
    mem->ptr_qp_in = qp_in;

    // reduce eq constr DOF
    d_ocp_qp_reduce_eq_dof(qp_in, mem->red_qp, opts->hpipm_red_opts, mem->hpipm_red_work);

    // condense gradient only
    d_cond_qp_cond_rhs(mem->red_qp, fcond_qp_in, opts->hpipm_cond_opts, mem->hpipm_cond_work);

    // stop timer
    mem->time_qp_xcond = acados_toc(&timer);

    return ACADOS_SUCCESS;
}



int ocp_qp_full_expansion(void *fcond_qp_out_, void *qp_out_, void *opts_, void *mem_, void *work)
{
    dense_qp_out *fcond_qp_out = fcond_qp_out_;
    ocp_qp_out *qp_out = qp_out_;
    ocp_qp_full_condensing_opts *opts = opts_;
    ocp_qp_full_condensing_memory *mem = mem_;

    acados_timer timer;

    // start timer
    acados_tic(&timer);

    // expand solution
    if (opts->expand_dual_sol == 0)
    {
        d_cond_qp_expand_primal_sol(mem->red_qp, fcond_qp_out, mem->red_sol, opts->hpipm_cond_opts, mem->hpipm_cond_work);
    }
    else
    {
        d_cond_qp_expand_sol(mem->red_qp, fcond_qp_out, mem->red_sol, opts->hpipm_cond_opts, mem->hpipm_cond_work);
    }

    // restore solution
    d_ocp_qp_restore_eq_dof(mem->ptr_qp_in, mem->red_sol, qp_out, opts->hpipm_red_opts, mem->hpipm_red_work);

//d_ocp_qp_sol_print(mem->red_sol->dim, mem->red_sol);
//d_ocp_qp_sol_print(qp_out->dim, qp_out);
//exit(1);
    // stop timer
    mem->time_qp_xcond += acados_toc(&timer);

    return ACADOS_SUCCESS;
}



void ocp_qp_full_condensing_config_initialize_default(void *config_)
{
    ocp_qp_xcond_config *config = config_;

    config->dims_calculate_size = &ocp_qp_full_condensing_dims_calculate_size;
    config->dims_assign = &ocp_qp_full_condensing_dims_assign;
    config->dims_set = &ocp_qp_full_condensing_dims_set;
    config->dims_get = &ocp_qp_full_condensing_dims_get;
    config->opts_calculate_size = &ocp_qp_full_condensing_opts_calculate_size;
    config->opts_assign = &ocp_qp_full_condensing_opts_assign;
    config->opts_initialize_default = &ocp_qp_full_condensing_opts_initialize_default;
    config->opts_update = &ocp_qp_full_condensing_opts_update;
    config->opts_set = &ocp_qp_full_condensing_opts_set;
    config->memory_calculate_size = &ocp_qp_full_condensing_memory_calculate_size;
    config->memory_assign = &ocp_qp_full_condensing_memory_assign;
    config->memory_get = &ocp_qp_full_condensing_memory_get;
    config->workspace_calculate_size = &ocp_qp_full_condensing_workspace_calculate_size;
    config->condensing = &ocp_qp_full_condensing;
    config->condensing_rhs = &ocp_qp_full_condensing_rhs;
    config->expansion = &ocp_qp_full_expansion;

    return;
}
