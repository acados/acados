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


// external
#include <stdlib.h>
#include <assert.h>
#include <string.h>
// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing.h"
#include "acados/utils/mem.h"
#include "hpipm/include/hpipm_d_ocp_qp_red.h"
// hpipm
#include "hpipm/include/hpipm_d_cond.h"
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_dim.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_part_cond.h"
#include "acados/utils/timing.h"



/************************************************
 * dims
 ************************************************/

acados_size_t ocp_qp_partial_condensing_dims_calculate_size(void *config_, int N)
{
    acados_size_t size = 0;

    size += sizeof(ocp_qp_partial_condensing_dims);

    // orig_dims, red_dims
    size += 2*sizeof(ocp_qp_dims);
    size += 2*d_ocp_qp_dim_memsize(N);

    // pcond_dims
    size += sizeof(ocp_qp_dims);
    size += d_ocp_qp_dim_memsize(N);  // worst-case size of new QP

    // block_size
    size += (N + 1) * sizeof(int);

    size += 2*8;

    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_qp_partial_condensing_dims_assign(void *config_, int N, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // dims
    ocp_qp_partial_condensing_dims *dims = (ocp_qp_partial_condensing_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_partial_condensing_dims);

    // orig_dims
    dims->orig_dims = (ocp_qp_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_dims);
    // red_dims
    dims->red_dims = (ocp_qp_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_dims);
    // pcond_dims
    dims->pcond_dims = (ocp_qp_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_dims);

    align_char_to(8, &c_ptr);

    // orig_dims
    d_ocp_qp_dim_create(N, dims->orig_dims, c_ptr);
    c_ptr += d_ocp_qp_dim_memsize(N);
    // red_dims
    d_ocp_qp_dim_create(N, dims->red_dims, c_ptr);
    c_ptr += d_ocp_qp_dim_memsize(N);
    // pcond_dims
    d_ocp_qp_dim_create(N, dims->pcond_dims, c_ptr);
    c_ptr += d_ocp_qp_dim_memsize(N);

    // block_size
    assign_and_advance_int(N + 1, &dims->block_size, &c_ptr);

    align_char_to(8, &c_ptr);

    assert((char *) raw_memory + ocp_qp_partial_condensing_dims_calculate_size(config_, N) >= c_ptr);

    return dims;
}



void ocp_qp_partial_condensing_dims_set(void *config_, void *dims_, int stage, const char *field, int *value)
{
    ocp_qp_partial_condensing_dims *dims = dims_;

    ocp_qp_dims_set(config_, dims->orig_dims, stage, field, value);

    // TODO later in mem do the pcond_dims

    return;
}



void ocp_qp_partial_condensing_dims_get(void *config_, void *dims_, const char *field, void* value)
{
    ocp_qp_partial_condensing_dims *dims = dims_;

    if(!strcmp(field, "xcond_dims"))
    {
        ocp_qp_dims **ptr = value;
        *ptr = dims->pcond_dims;
    }
    else
    {
        printf("\nerror: ocp_qp_partial_condensing_dims_get: field %s not available\n", field);
        exit(1);
    }

    return;

}



/************************************************
 * opts
 ************************************************/

acados_size_t ocp_qp_partial_condensing_opts_calculate_size(void *dims_)
{
    ocp_qp_partial_condensing_dims *dims = dims_;

    size_t N = dims->orig_dims->N;

    // populate dimensions of reduced qp
    d_ocp_qp_dim_reduce_eq_dof(dims->orig_dims, dims->red_dims);

    // (temporarily) populate dimensions of new ocp_qp based on N2==N
    // TODO: needed?
    size_t N2 = N;
    dims->pcond_dims->N = N2;
    d_part_cond_qp_compute_block_size(dims->red_dims->N, N2, dims->block_size);
    d_part_cond_qp_compute_dim(dims->red_dims, dims->block_size, dims->pcond_dims);

    acados_size_t size = 0;

    size += sizeof(ocp_qp_partial_condensing_opts);

    // block_size
    size += (N + 1) * sizeof(int);

    // hpipm_pcond_opts
    size += sizeof(struct d_part_cond_qp_arg);
    size += d_part_cond_qp_arg_memsize(N);  // worst case size of new QP
    // hpipm_red_opts
    size += sizeof(struct d_ocp_qp_reduce_eq_dof_arg);
    size += d_ocp_qp_reduce_eq_dof_arg_memsize();

    size += 3*8;
    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_qp_partial_condensing_opts_assign(void *dims_, void *raw_memory)
{
    ocp_qp_partial_condensing_dims *dims = dims_;

    int N = dims->orig_dims->N;

    char *c_ptr = (char *) raw_memory;

    // opts
    ocp_qp_partial_condensing_opts *opts = (ocp_qp_partial_condensing_opts *) c_ptr;
    c_ptr += sizeof(ocp_qp_partial_condensing_opts);

    align_char_to(8, &c_ptr);

    // hpipm_pcond_opts
    opts->hpipm_pcond_opts = (struct d_part_cond_qp_arg *) c_ptr;
    c_ptr += sizeof(struct d_part_cond_qp_arg);
    // hpipm_red_opts
    opts->hpipm_red_opts = (struct d_ocp_qp_reduce_eq_dof_arg *) c_ptr;
    c_ptr += sizeof(struct d_ocp_qp_reduce_eq_dof_arg);

    align_char_to(8, &c_ptr);

    // hpipm_pcond_opts
    d_part_cond_qp_arg_create(N, opts->hpipm_pcond_opts, c_ptr);
    c_ptr += opts->hpipm_pcond_opts->memsize;
    // hpipm_red_opts
    d_ocp_qp_reduce_eq_dof_arg_create(opts->hpipm_red_opts, c_ptr);
    c_ptr += opts->hpipm_red_opts->memsize;

    // block_size
    assign_and_advance_int(N + 1, &opts->block_size, &c_ptr);
    align_char_to(8, &c_ptr);
    opts->block_size_was_set = false;

    assert((char *) raw_memory + ocp_qp_partial_condensing_opts_calculate_size(dims) >= c_ptr);

    return opts;
}



void ocp_qp_partial_condensing_opts_initialize_default(void *dims_, void *opts_)
{
    ocp_qp_partial_condensing_dims *dims = dims_;
    ocp_qp_partial_condensing_opts *opts = opts_;

    int N = dims->orig_dims->N;

    opts->N2 = N;  // no partial condensing by default
    opts->N2_bkp = opts->N2;

    dims->pcond_dims->N = opts->N2;
    opts->hpipm_pcond_opts->N2 = opts->N2;
    // hpipm_pcond_opts
    d_part_cond_qp_arg_set_default(opts->hpipm_pcond_opts);
    // hpipm_red_opts
    int tmp_i1 = 1;
    d_ocp_qp_reduce_eq_dof_arg_set_default(opts->hpipm_red_opts);
    d_ocp_qp_reduce_eq_dof_arg_set_alias_unchanged(&tmp_i1, opts->hpipm_red_opts);

    opts->mem_qp_in = 1;

    return;
}



void ocp_qp_partial_condensing_opts_update(void *dims_, void *opts_)
{
    ocp_qp_partial_condensing_dims *dims = dims_;
    ocp_qp_partial_condensing_opts *opts = opts_;

    dims->pcond_dims->N = opts->N2;
    opts->N2_bkp = opts->N2;
    // hpipm_pcond_opts
    opts->hpipm_pcond_opts->N2 = opts->N2;
    d_part_cond_qp_arg_set_default(opts->hpipm_pcond_opts);
    d_part_cond_qp_arg_set_ric_alg(opts->ric_alg, opts->hpipm_pcond_opts);

    return;
}



void ocp_qp_partial_condensing_opts_set(void *opts_, const char *field, void* value)
{
    ocp_qp_partial_condensing_opts *opts = opts_;

    if(!strcmp(field, "N"))
    {
        int *tmp_ptr = value;
        opts->N2 = *tmp_ptr;
    }
    else if(!strcmp(field, "N_bkp"))
    {
        int *tmp_ptr = value;
        opts->N2_bkp = *tmp_ptr;
    }
    else if(!strcmp(field, "ric_alg"))
    {
        int *tmp_ptr = value;
        opts->ric_alg = *tmp_ptr;
    }
    else if(!strcmp(field, "block_size"))
    {
        int *tmp_ptr = value;
        for (int i = 0; i < opts->N2+1; i++)
        {
            opts->block_size[i] = tmp_ptr[i];
        }
        opts->block_size_was_set = true;
    }
    // TODO dual_sol ???
    else
    {
        printf("\nerror: field %s not available in ocp_qp_partial_condensing_opts_set\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_qp_partial_condensing_memory_calculate_size(void *dims_, void *opts_)
{
    ocp_qp_partial_condensing_opts *opts = opts_;
    ocp_qp_partial_condensing_dims *dims = dims_;

    // populate dimensions of reduced qp
    d_ocp_qp_dim_reduce_eq_dof(dims->orig_dims, dims->red_dims);

    // populate dimensions from options: N2, block_size
    dims->pcond_dims->N = opts->N2;
    if (opts->block_size_was_set)
    {
        for (int i=0; i<opts->N2+1; i++)
        {
            dims->block_size[i] = opts->block_size[i];
        }
        // check consistency N and sum(block_size)
        int tmp_sum = 0;
        for (int i = 0; i<opts->N2+1; i++)
        {
            tmp_sum += dims->block_size[i];
        }
        if (tmp_sum != dims->orig_dims->N)
        {
            printf("partial condensing: sum of block_size should match N, got %d != N = %d\n", tmp_sum, dims->orig_dims->N);
            exit(1);
        }
    }
    else
    {
        d_part_cond_qp_compute_block_size(dims->red_dims->N, opts->N2, dims->block_size);
    }

    // printf("using block_size:\n");
    // tmp_sum = 0;
    // for (int i = 0; i<opts->N2+1; i++)
    // {
    //     tmp_sum += dims->block_size[i];
    //     printf("%d\t", dims->block_size[i]);
    // }
    // printf("\ntotal %d\n", tmp_sum);

    d_part_cond_qp_compute_dim(dims->red_dims, dims->block_size, dims->pcond_dims);

    // printf("\n\ndimensions of partially condensed QP:\n");
    // print_ocp_qp_dims(dims->pcond_dims);

    // printf("\n\ndimensions of original QP:\n");
    // print_ocp_qp_dims(dims->orig_dims);

    acados_size_t size = 0;

    size += sizeof(ocp_qp_partial_condensing_memory);

    size += ocp_qp_in_calculate_size(dims->pcond_dims);
    size += ocp_qp_out_calculate_size(dims->pcond_dims);
    size += ocp_qp_seed_calculate_size(dims->pcond_dims);

    size += ocp_qp_in_calculate_size(dims->red_dims);
    size += ocp_qp_out_calculate_size(dims->red_dims);
    size += ocp_qp_seed_calculate_size(dims->red_dims);

    // hpipm_pcond_work
    size += sizeof(struct d_part_cond_qp_ws);
    size += d_part_cond_qp_ws_memsize(dims->red_dims, dims->block_size, dims->pcond_dims, opts->hpipm_pcond_opts);

    size += sizeof(struct d_ocp_qp_reduce_eq_dof_ws);
    size += d_ocp_qp_reduce_eq_dof_ws_memsize(dims->orig_dims);

    size += 2*8;
    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_qp_partial_condensing_memory_assign(void *dims_, void *opts_, void *raw_memory)
{
    ocp_qp_partial_condensing_dims *dims = dims_;
    ocp_qp_partial_condensing_opts *opts = opts_;

    char *c_ptr = (char *) raw_memory;

    // initial alignment
    align_char_to(8, &c_ptr);

    ocp_qp_partial_condensing_memory *mem = (ocp_qp_partial_condensing_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_partial_condensing_memory);

    align_char_to(8, &c_ptr);

    // hpipm_pcond_work struct
    mem->hpipm_pcond_work = (struct d_part_cond_qp_ws *) c_ptr;
    c_ptr += sizeof(struct d_part_cond_qp_ws);
    // hpipm_red_work struct
    mem->hpipm_red_work = (struct d_ocp_qp_reduce_eq_dof_ws *) c_ptr;
    c_ptr += sizeof(struct d_ocp_qp_reduce_eq_dof_ws);
    align_char_to(8, &c_ptr);

    // hpipm_pcond_work
    d_part_cond_qp_ws_create(dims->red_dims, dims->block_size, dims->pcond_dims,
                             opts->hpipm_pcond_opts, mem->hpipm_pcond_work, c_ptr);
    c_ptr += mem->hpipm_pcond_work->memsize;
    // hpipm_red_work
    d_ocp_qp_reduce_eq_dof_ws_create(dims->orig_dims, mem->hpipm_red_work, c_ptr);
    c_ptr += mem->hpipm_red_work->memsize;

    mem->pcond_qp_in = ocp_qp_in_assign(dims->pcond_dims, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(dims->pcond_dims);

    mem->pcond_qp_out = ocp_qp_out_assign(dims->pcond_dims, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(dims->pcond_dims);

    mem->pcond_qp_seed = ocp_qp_seed_assign(dims->pcond_dims, c_ptr);
    c_ptr += ocp_qp_seed_calculate_size(dims->pcond_dims);

    mem->red_qp = ocp_qp_in_assign(dims->red_dims, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(dims->red_dims);

    mem->red_sol = ocp_qp_out_assign(dims->red_dims, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(dims->red_dims);

    mem->red_seed = ocp_qp_seed_assign(dims->red_dims, c_ptr);
    c_ptr += ocp_qp_seed_calculate_size(dims->red_dims);

    mem->qp_out_info = (qp_info *) mem->pcond_qp_out->misc;

    mem->dims = dims;

    assert((char *) raw_memory + ocp_qp_partial_condensing_memory_calculate_size(dims, opts) >= c_ptr);

    return mem;
}



void ocp_qp_partial_condensing_memory_get(void *config_, void *mem_, const char *field, void* value)
{
    ocp_qp_partial_condensing_memory *mem = mem_;

    if(!strcmp(field, "xcond_qp_in"))
    {
        ocp_qp_in **ptr = value;
        *ptr = mem->pcond_qp_in;
    }
    else if(!strcmp(field, "xcond_qp_out"))
    {
        ocp_qp_out **ptr = value;
        *ptr = mem->pcond_qp_out;
    }
    else if(!strcmp(field, "xcond_seed"))
    {
        ocp_qp_seed **ptr = value;
        *ptr = mem->pcond_qp_seed;
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
        printf("\nerror: ocp_qp_partial_condensing_memory_get: field %s not available\n", field);
        exit(1);
    }

    return;

}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_qp_partial_condensing_workspace_calculate_size(void *dims_, void *opts_)
{
    return 0;
}



/************************************************
 * functions
 ************************************************/

int ocp_qp_partial_condensing(void *qp_in_, void *pcond_qp_in_, void *opts_, void *mem_, void *work)
{
    ocp_qp_in *qp_in = qp_in_;
    ocp_qp_in *pcond_qp_in = pcond_qp_in_;
    ocp_qp_partial_condensing_opts *opts = opts_;
    ocp_qp_partial_condensing_memory *mem = mem_;

    acados_timer timer;

    assert(opts->N2 == opts->N2_bkp);

    // save pointers to ocp_qp_in in memory (needed for expansion)
    mem->ptr_qp_in = qp_in;
    mem->ptr_pcond_qp_in = pcond_qp_in;

    // start timer
    acados_tic(&timer);

    // reduce eq constr DOF
    d_ocp_qp_reduce_eq_dof(qp_in, mem->red_qp, opts->hpipm_red_opts, mem->hpipm_red_work);

    // d_ocp_qp_codegen("orig_ocp_qp.c", "w", qp_in->dim, qp_in);
    // d_ocp_qp_codegen("pcond_ocp_qp.c", "w", pcond_qp_in->dim, pcond_qp_in);
    // d_ocp_qp_print(qp_in->dim, qp_in);
    // d_ocp_qp_print(pcond_qp_in->dim, pcond_qp_in);

    // convert to partially condensed qp structure
    d_part_cond_qp_cond(mem->red_qp, pcond_qp_in, opts->hpipm_pcond_opts, mem->hpipm_pcond_work);

    // stop timer
    mem->time_qp_xcond = acados_toc(&timer);

    return ACADOS_SUCCESS;
}


int ocp_qp_partial_condensing_condense_qp_out(void *qp_in_, void *pcond_qp_in_, void *qp_out_, void *pcond_qp_out_, void *opts_, void *mem_, void *work)
{
    ocp_qp_in *qp_in = qp_in_;
    ocp_qp_in *pcond_qp_in = pcond_qp_in_;
    ocp_qp_out *qp_out = qp_out_;
    ocp_qp_out *pcond_qp_out = pcond_qp_out_;
    ocp_qp_partial_condensing_opts *opts = opts_;
    ocp_qp_partial_condensing_memory *mem = mem_;

    d_ocp_qp_reduce_eq_dof_sol(qp_in, qp_out, mem->red_sol, opts->hpipm_red_opts, mem->hpipm_red_work);
    d_part_cond_qp_cond_sol(mem->red_qp, pcond_qp_in, mem->red_sol, pcond_qp_out, opts->hpipm_pcond_opts, mem->hpipm_pcond_work);

    return ACADOS_SUCCESS;
}


int ocp_qp_partial_condensing_condense_lhs(void *qp_in_, void *pcond_qp_in_, void *opts_, void *mem_, void *work)
{
    ocp_qp_in *qp_in = qp_in_;
    ocp_qp_in *pcond_qp_in = pcond_qp_in_;
    ocp_qp_partial_condensing_opts *opts = opts_;
    ocp_qp_partial_condensing_memory *mem = mem_;

    acados_timer timer;

    assert(opts->N2 == opts->N2_bkp);

    // save pointers to ocp_qp_in in memory (needed for expansion)
    mem->ptr_qp_in = qp_in;
    mem->ptr_pcond_qp_in = pcond_qp_in;

    acados_tic(&timer);

    d_ocp_qp_reduce_eq_dof_lhs(qp_in, mem->red_qp, opts->hpipm_red_opts, mem->hpipm_red_work);
    d_part_cond_qp_cond_lhs(mem->red_qp, pcond_qp_in, opts->hpipm_pcond_opts, mem->hpipm_pcond_work);

    mem->time_qp_xcond = acados_toc(&timer);

    return ACADOS_SUCCESS;
}



int ocp_qp_partial_condensing_condense_rhs(void *qp_in_, void *pcond_qp_in_, void *opts_, void *mem_, void *work)
{
    ocp_qp_in *qp_in = qp_in_;
    ocp_qp_in *pcond_qp_in = pcond_qp_in_;
    ocp_qp_partial_condensing_opts *opts = opts_;
    ocp_qp_partial_condensing_memory *mem = mem_;

    assert(opts->N2 == opts->N2_bkp);

    acados_timer timer;

    // start timer
    acados_tic(&timer);

    // save pointers to ocp_qp_in in memory (needed for expansion)
    mem->ptr_qp_in = qp_in;
    mem->ptr_pcond_qp_in = pcond_qp_in;

    // reduce eq constr DOF
    d_ocp_qp_reduce_eq_dof_rhs(qp_in, mem->red_qp, opts->hpipm_red_opts, mem->hpipm_red_work);

    // convert to partially condensed qp structure
    d_part_cond_qp_cond_rhs(mem->red_qp, pcond_qp_in, opts->hpipm_pcond_opts, mem->hpipm_pcond_work);

    // stop timer
    mem->time_qp_xcond += acados_toc(&timer);

    return ACADOS_SUCCESS;
}



int ocp_qp_partial_condensing_condense_rhs_seed(void *qp_in_, void *qp_seed, void *pcond_seed, void *opts_, void *mem_, void *work)
{
    ocp_qp_in *qp_in = qp_in_;
    ocp_qp_partial_condensing_opts *opts = opts_;
    ocp_qp_partial_condensing_memory *mem = mem_;

    assert(opts->N2 == opts->N2_bkp);

    acados_timer timer;

    // start timer
    acados_tic(&timer);

    // save pointers to ocp_qp_in in memory (needed for expansion)
    mem->ptr_qp_in = qp_in;
    mem->ptr_qp_seed = qp_seed;

    // reduce eq constr DOF: residual
    d_ocp_qp_reduce_eq_dof_seed(qp_in, qp_seed, mem->red_seed, opts->hpipm_red_opts, mem->hpipm_red_work);

    // convert to partially condensed qp structure
    d_part_cond_qp_cond_seed(mem->red_qp, mem->red_seed, pcond_seed, opts->hpipm_pcond_opts, mem->hpipm_pcond_work);

    // stop timer
    mem->time_qp_xcond += acados_toc(&timer);

    return ACADOS_SUCCESS;
}


int ocp_qp_partial_expansion(void *pcond_qp_out_, void *qp_out_, void *opts_, void *mem_, void *work)
{
    ocp_qp_out *pcond_qp_out = pcond_qp_out_;
    ocp_qp_out *qp_out = qp_out_;
    ocp_qp_partial_condensing_opts *opts = opts_;
    ocp_qp_partial_condensing_memory *mem = mem_;

    assert(opts->N2 == opts->N2_bkp);

    acados_timer timer;

    // start timer
    acados_tic(&timer);

    // expand solution
    // TODO only if N2<N
    d_part_cond_qp_expand_sol(mem->red_qp, pcond_qp_out, mem->red_sol, opts->hpipm_pcond_opts, mem->hpipm_pcond_work);

    // restore solution
    d_ocp_qp_restore_eq_dof(mem->ptr_qp_in, mem->red_sol, qp_out, opts->hpipm_red_opts, mem->hpipm_red_work);

    // stop timer
    mem->time_qp_xcond += acados_toc(&timer);

    return ACADOS_SUCCESS;
}

int ocp_qp_partial_condensing_expand_sol_seed(void *pcond_qp_out_, void *qp_out_, void *opts_, void *mem_, void *work)
{
    ocp_qp_out *pcond_qp_out = pcond_qp_out_;
    ocp_qp_out *qp_out = qp_out_;
    ocp_qp_partial_condensing_opts *opts = opts_;
    ocp_qp_partial_condensing_memory *mem = mem_;

    assert(opts->N2 == opts->N2_bkp);

    acados_timer timer;

    // start timer
    acados_tic(&timer);

    // expand solution
    d_part_cond_qp_expand_sol_seed(mem->red_qp, mem->red_seed, pcond_qp_out, mem->red_sol, opts->hpipm_pcond_opts, mem->hpipm_pcond_work);

    // restore solution
    d_ocp_qp_restore_eq_dof_seed(mem->ptr_qp_in, mem->ptr_qp_seed, mem->red_sol, qp_out, opts->hpipm_red_opts, mem->hpipm_red_work);

    // stop timer
    mem->time_qp_xcond += acados_toc(&timer);

    return ACADOS_SUCCESS;
}



void ocp_qp_partial_condensing_config_initialize_default(void *config_)
{
    ocp_qp_xcond_config *config = config_;

    config->dims_calculate_size = &ocp_qp_partial_condensing_dims_calculate_size;
    config->dims_assign = &ocp_qp_partial_condensing_dims_assign;
    config->dims_set = &ocp_qp_partial_condensing_dims_set;
    config->dims_get = &ocp_qp_partial_condensing_dims_get;
    config->opts_calculate_size = &ocp_qp_partial_condensing_opts_calculate_size;
    config->opts_assign = &ocp_qp_partial_condensing_opts_assign;
    config->opts_initialize_default = &ocp_qp_partial_condensing_opts_initialize_default;
    config->opts_update = &ocp_qp_partial_condensing_opts_update;
    config->opts_set = &ocp_qp_partial_condensing_opts_set;
    config->memory_calculate_size = &ocp_qp_partial_condensing_memory_calculate_size;
    config->memory_assign = &ocp_qp_partial_condensing_memory_assign;
    config->memory_get = &ocp_qp_partial_condensing_memory_get;
    config->workspace_calculate_size = &ocp_qp_partial_condensing_workspace_calculate_size;
    config->condensing = &ocp_qp_partial_condensing;
    config->condense_lhs = &ocp_qp_partial_condensing_condense_lhs;
    config->condense_rhs = &ocp_qp_partial_condensing_condense_rhs;
    config->condense_rhs_seed = &ocp_qp_partial_condensing_condense_rhs_seed;
    config->condense_qp_out = &ocp_qp_partial_condensing_condense_qp_out;
    config->expansion = &ocp_qp_partial_expansion;
    config->expand_sol_seed = &ocp_qp_partial_condensing_expand_sol_seed;

    return;
}
