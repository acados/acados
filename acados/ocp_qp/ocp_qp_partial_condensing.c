/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

// external
#include <assert.h>
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

/************************************************
 * opts
 ************************************************/

int ocp_qp_partial_condensing_opts_calculate_size(ocp_qp_dims *dims)
{
    int N = dims->N;

    int size = 0;

    size += sizeof(ocp_qp_partial_condensing_opts);

    // hpipm opts
    size += sizeof(struct d_cond_qp_ocp2ocp_arg);
    size += d_memsize_cond_qp_ocp2ocp_arg(N);  // worst case size of new QP
                                               //
    size += sizeof(ocp_qp_dims);
    size += d_memsize_ocp_qp_dim(N);  // worst-case size of new QP
    size += (N + 1) * sizeof(int);    // block size

    size += 1 * 8;
    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_qp_partial_condensing_opts_assign(ocp_qp_dims *dims, void *raw_memory)
{
    int N = dims->N;

    char *c_ptr = (char *) raw_memory;

    // opts
    ocp_qp_partial_condensing_opts *opts = (ocp_qp_partial_condensing_opts *) c_ptr;
    c_ptr += sizeof(ocp_qp_partial_condensing_opts);

    // pcond_dims
    opts->pcond_dims = (ocp_qp_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_dims);
    // hpipm_opts
    opts->hpipm_opts = (struct d_cond_qp_ocp2ocp_arg *) c_ptr;
    c_ptr += sizeof(struct d_cond_qp_ocp2ocp_arg);

    // block size
    assign_and_advance_int(N + 1, &opts->block_size, &c_ptr);

    align_char_to(8, &c_ptr);

    // pcond_dims
    d_create_ocp_qp_dim(N, opts->pcond_dims, c_ptr);
    c_ptr += d_memsize_ocp_qp_dim(dims->N);
    // hpipm_opts
    d_create_cond_qp_ocp2ocp_arg(N, opts->hpipm_opts, c_ptr);
    c_ptr += opts->hpipm_opts->memsize;

    assert((char *) raw_memory + ocp_qp_partial_condensing_opts_calculate_size(dims) >= c_ptr);

    return opts;
}



void ocp_qp_partial_condensing_opts_initialize_default(ocp_qp_dims *dims, void *opts_)
{
    ocp_qp_partial_condensing_opts *opts = opts_;

    int N = dims->N;

    opts->N2 = N;  // no partial condensing by default
    opts->N2_bkp = opts->N2;

    opts->pcond_dims->N = opts->N2;
    // hpipm_opts
    d_set_default_cond_qp_ocp2ocp_arg(opts->N2, opts->hpipm_opts);
}



void ocp_qp_partial_condensing_opts_update(ocp_qp_dims *dims, void *opts_)
{
    ocp_qp_partial_condensing_opts *opts = opts_;

    opts->pcond_dims->N = opts->N2;
    opts->N2_bkp = opts->N2;
    // hpipm_opts
    d_set_default_cond_qp_ocp2ocp_arg(opts->N2, opts->hpipm_opts);
}

/************************************************
 * memory
 ************************************************/

int ocp_qp_partial_condensing_memory_calculate_size(ocp_qp_dims *dims, void *opts_)
{
    ocp_qp_partial_condensing_opts *opts = opts_;

    int size = 0;

    // populate dimensions of new ocp_qp based on N2
    opts->pcond_dims->N = opts->N2;
    // TODO(all): user-defined block size
    d_compute_block_size_cond_qp_ocp2ocp(dims->N, opts->N2, opts->block_size);

    d_compute_qp_dim_ocp2ocp(dims, opts->block_size, opts->pcond_dims);

    size += sizeof(ocp_qp_partial_condensing_memory);
    // hpipm workspace
    size += sizeof(struct d_cond_qp_ocp2ocp_workspace);
    size += d_memsize_cond_qp_ocp2ocp(dims, opts->block_size, opts->pcond_dims, opts->hpipm_opts);

    size += 1 * 8;
    return size;
}

void *ocp_qp_partial_condensing_memory_assign(ocp_qp_dims *dims, void *opts_, void *raw_memory)
{
    ocp_qp_partial_condensing_opts *opts = opts_;

    char *c_ptr = (char *) raw_memory;

    ocp_qp_partial_condensing_memory *mem = (ocp_qp_partial_condensing_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_partial_condensing_memory);
    // hpipm_workspace
    mem->hpipm_workspace = (struct d_cond_qp_ocp2ocp_workspace *) c_ptr;
    c_ptr += sizeof(struct d_cond_qp_ocp2ocp_workspace);

    // hpipm workspace structure
    align_char_to(8, &c_ptr);
    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    // hpipm_workspace
    d_create_cond_qp_ocp2ocp(dims, opts->block_size, opts->pcond_dims, opts->hpipm_opts,
                             mem->hpipm_workspace, c_ptr);
    c_ptr += mem->hpipm_workspace->memsize;

    mem->qp_in = NULL;  // initialized when partial condensing routine is called

    assert((char *) raw_memory + ocp_qp_partial_condensing_memory_calculate_size(dims, opts) >=
           c_ptr);

    return mem;
}

/************************************************
 * memory
 ************************************************/

int ocp_qp_partial_condensing_workspace_calculate_size(ocp_qp_dims *dims, void *opts_) { return 0; }
/************************************************
 * functions
 ************************************************/

void ocp_qp_partial_condensing(ocp_qp_in *in, ocp_qp_in *out, ocp_qp_partial_condensing_opts *opts,
                               ocp_qp_partial_condensing_memory *mem, void *work)
{
    assert(opts->N2 == opts->N2_bkp);

    // save pointers to ocp_qp_in in memory (needed for expansion)
    mem->qp_in = in;
    mem->pcond_qp_in = out;

    // convert to partially condensed qp structure
    d_cond_qp_ocp2ocp(in, out, opts->hpipm_opts, mem->hpipm_workspace);
}

void ocp_qp_partial_expansion(ocp_qp_out *in, ocp_qp_out *out, ocp_qp_partial_condensing_opts *opts,
                              ocp_qp_partial_condensing_memory *mem, void *work)
{
    assert(opts->N2 == opts->N2_bkp);

    d_expand_sol_ocp2ocp(mem->qp_in, mem->pcond_qp_in, in, out, opts->hpipm_opts,
                         mem->hpipm_workspace);
}

void ocp_qp_partial_condensing_config_initialize_default(void *config_)
{
    ocp_qp_condensing_config *config = config_;

    config->opts_calculate_size = &ocp_qp_partial_condensing_opts_calculate_size;
    config->opts_assign = &ocp_qp_partial_condensing_opts_assign;
    config->opts_initialize_default = &ocp_qp_partial_condensing_opts_initialize_default;
    config->opts_update = &ocp_qp_partial_condensing_opts_update;
    config->memory_calculate_size = &ocp_qp_partial_condensing_memory_calculate_size;
    config->memory_assign = &ocp_qp_partial_condensing_memory_assign;
    config->workspace_calculate_size = &ocp_qp_partial_condensing_workspace_calculate_size;
    // TODO(dimitris): either do casting as below or void in defs, as above (also pass config or
    // not?)
    config->condensing =
        (int (*)(void *, void *, void *, void *, void *)) & ocp_qp_partial_condensing;
    config->expansion =
        (int (*)(void *, void *, void *, void *, void *)) & ocp_qp_partial_expansion;

    return;
}
