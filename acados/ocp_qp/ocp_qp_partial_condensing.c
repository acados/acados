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
// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing.h"
#include "acados/utils/mem.h"
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

int ocp_qp_partial_condensing_dims_calculate_size(void *config_, int N)
{
    int size = 0;

    size += sizeof(ocp_qp_partial_condensing_dims);

	// orig_dims
    size += sizeof(ocp_qp_dims);
    size += d_ocp_qp_dim_memsize(N);

	// pcond_dims
    size += sizeof(ocp_qp_dims);
    size += d_ocp_qp_dim_memsize(N);  // worst-case size of new QP

	// block size
    size += (N + 1) * sizeof(int);

    size += 1 * 8;

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
    // pcond_dims
    dims->pcond_dims = (ocp_qp_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_dims);

    // block size
    assign_and_advance_int(N + 1, &dims->block_size, &c_ptr);

    align_char_to(8, &c_ptr);

    // orig_dims
    d_ocp_qp_dim_create(N, dims->orig_dims, c_ptr);
    c_ptr += d_ocp_qp_dim_memsize(N);
    // pcond_dims
    d_ocp_qp_dim_create(N, dims->pcond_dims, c_ptr);
    c_ptr += d_ocp_qp_dim_memsize(N);

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

int ocp_qp_partial_condensing_opts_calculate_size(void *dims_)
{
	ocp_qp_partial_condensing_dims *dims = dims_;

    int N = dims->orig_dims->N;

    // (temporarely) populate dimensions of new ocp_qp based on N2==N
	int N2 = N;
    dims->pcond_dims->N = N2;
    // TODO(all): user-defined block size
    d_part_cond_qp_compute_block_size(dims->orig_dims->N, N2, dims->block_size);
    d_part_cond_qp_compute_dim(dims->orig_dims, dims->block_size, dims->pcond_dims);

    int size = 0;

    size += sizeof(ocp_qp_partial_condensing_opts);

    // hpipm opts
    size += sizeof(struct d_part_cond_qp_arg);
    size += d_part_cond_qp_arg_memsize(N);  // worst case size of new QP
                                               //
	// pcond_dims
//    size += sizeof(ocp_qp_dims);
//    size += d_ocp_qp_dim_memsize(N);  // worst-case size of new QP

	// block size
//    size += (N + 1) * sizeof(int);

    size += 1 * 8;
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

    // pcond_dims
//    opts->pcond_dims = (ocp_qp_dims *) c_ptr;
//    c_ptr += sizeof(ocp_qp_dims);

    // hpipm_opts
    opts->hpipm_opts = (struct d_part_cond_qp_arg *) c_ptr;
    c_ptr += sizeof(struct d_part_cond_qp_arg);

    // block size
//    assign_and_advance_int(N + 1, &opts->block_size, &c_ptr);

    align_char_to(8, &c_ptr);

    // pcond_dims
//    d_ocp_qp_dim_create(N, opts->pcond_dims, c_ptr);
//    c_ptr += d_ocp_qp_dim_memsize(dims->N);
    // hpipm_opts
    d_part_cond_qp_arg_create(N, opts->hpipm_opts, c_ptr);
    c_ptr += opts->hpipm_opts->memsize;

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
    // hpipm_opts
    d_part_cond_qp_arg_set_default(opts->N2, opts->hpipm_opts);

	opts->mem_qp_in = 1;

	return;
}



void ocp_qp_partial_condensing_opts_update(void *dims_, void *opts_)
{
	ocp_qp_partial_condensing_dims *dims = dims_;
    ocp_qp_partial_condensing_opts *opts = opts_;

    dims->pcond_dims->N = opts->N2;
    opts->N2_bkp = opts->N2;
    // hpipm_opts
    d_part_cond_qp_arg_set_default(opts->N2, opts->hpipm_opts);
	d_part_cond_qp_arg_set_ric_alg(opts->ric_alg, opts->N2, opts->hpipm_opts);

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

int ocp_qp_partial_condensing_memory_calculate_size(void *dims_, void *opts_)
{
    ocp_qp_partial_condensing_opts *opts = opts_;
	ocp_qp_partial_condensing_dims *dims = dims_;

    // populate dimensions of new ocp_qp based on actual N2
    dims->pcond_dims->N = opts->N2;
    // TODO(all): user-defined block size
    d_part_cond_qp_compute_block_size(dims->orig_dims->N, opts->N2, dims->block_size);
    d_part_cond_qp_compute_dim(dims->orig_dims, dims->block_size, dims->pcond_dims);

    int size = 0;

    size += sizeof(ocp_qp_partial_condensing_memory);

	size += ocp_qp_in_calculate_size(dims->pcond_dims);

	size += ocp_qp_out_calculate_size(dims->pcond_dims);

    // hpipm workspace
    size += sizeof(struct d_part_cond_qp_ws);
    size += d_part_cond_qp_ws_memsize(dims->orig_dims, dims->block_size, dims->pcond_dims, opts->hpipm_opts);

    size += 2 * 8;

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

    // hpipm_workspace struct
    mem->hpipm_workspace = (struct d_part_cond_qp_ws *) c_ptr;
    c_ptr += sizeof(struct d_part_cond_qp_ws);

    // hpipm_workspace
    d_part_cond_qp_ws_create(dims->orig_dims, dims->block_size, dims->pcond_dims, opts->hpipm_opts, mem->hpipm_workspace, c_ptr);
    c_ptr += mem->hpipm_workspace->memsize;

	mem->pcond_qp_in = ocp_qp_in_assign(dims->pcond_dims, c_ptr);
	c_ptr += ocp_qp_in_calculate_size(dims->pcond_dims);

	mem->pcond_qp_out = ocp_qp_out_assign(dims->pcond_dims, c_ptr);
	c_ptr += ocp_qp_out_calculate_size(dims->pcond_dims);

	mem->qp_out_info = (qp_info *) mem->pcond_qp_out->misc;

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

int ocp_qp_partial_condensing_workspace_calculate_size(void *dims_, void *opts_)
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

    // convert to partially condensed qp structure
	// TODO only if N2<N
    d_part_cond_qp_cond(qp_in, pcond_qp_in, opts->hpipm_opts, mem->hpipm_workspace);

	// stop timer
    mem->time_qp_xcond = acados_toc(&timer);

	return ACADOS_SUCCESS;
}



int ocp_qp_partial_condensing_rhs(void *qp_in_, void *pcond_qp_in_, void *opts_, void *mem_, void *work)
{
	ocp_qp_in *qp_in = qp_in_;
	ocp_qp_in *pcond_qp_in = pcond_qp_in_;
	ocp_qp_partial_condensing_opts *opts = opts_;
	ocp_qp_partial_condensing_memory *mem = mem_;

    assert(opts->N2 == opts->N2_bkp);

    // save pointers to ocp_qp_in in memory (needed for expansion)
    mem->ptr_qp_in = qp_in;
    mem->ptr_pcond_qp_in = pcond_qp_in;

    // convert to partially condensed qp structure
	// TODO only if N2<N
    d_part_cond_qp_cond_rhs(qp_in, pcond_qp_in, opts->hpipm_opts, mem->hpipm_workspace);

	return ACADOS_SUCCESS;
}



int ocp_qp_partial_expansion(void *pcond_qp_out_, void *qp_out_, void *opts_, void *mem_, void *work)
{
	ocp_qp_out *pcond_qp_out = pcond_qp_out_;
	ocp_qp_out *qp_out = qp_out_;
	ocp_qp_partial_condensing_opts *opts = opts_;
	ocp_qp_partial_condensing_memory *mem = mem_;

    assert(opts->N2 == opts->N2_bkp);

	// TODO only if N2<N
    d_part_cond_qp_expand_sol(mem->ptr_qp_in, mem->ptr_pcond_qp_in, pcond_qp_out, qp_out, opts->hpipm_opts, mem->hpipm_workspace);

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
    config->condensing_rhs = &ocp_qp_partial_condensing_rhs;
    config->expansion = &ocp_qp_partial_expansion;

    return;
}
