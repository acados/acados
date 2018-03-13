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
#include "acados/ocp_qp/ocp_qp_partial_condensing.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/mem.h"
// hpipm
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp_dim.h"
#include "hpipm/include/hpipm_d_cond.h"
#include "hpipm/include/hpipm_d_part_cond.h"



int ocp_qp_partial_condensing_opts_calculate_size(ocp_qp_dims *dims)
{
    int size = 0;
    size += sizeof(ocp_qp_partial_condensing_opts);
    size += sizeof(ocp_qp_dims);
    size += d_memsize_ocp_qp_dim(dims->N);  // worst-case size of new QP
    return size;
}



ocp_qp_partial_condensing_opts *ocp_qp_partial_condensing_opts_assign(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_partial_condensing_opts *opts = (ocp_qp_partial_condensing_opts *)c_ptr;
    c_ptr += sizeof(ocp_qp_partial_condensing_opts);

    opts->pcond_dims = (ocp_qp_dims *)c_ptr;
    c_ptr += sizeof(ocp_qp_dims);

    assert((size_t)c_ptr % 8 == 0 && "double not 8-byte aligned!");

    d_create_ocp_qp_dim(dims->N, opts->pcond_dims, c_ptr);
    c_ptr += d_memsize_ocp_qp_dim(dims->N);

    assert((char*)raw_memory + ocp_qp_partial_condensing_opts_calculate_size(dims) == c_ptr);

    return opts;
}



void ocp_qp_partial_condensing_opts_initialize_default(ocp_qp_dims *dims, ocp_qp_partial_condensing_opts *opts)
{
    opts->N2 = dims->N;  // no partial condensing by default
    opts->pcond_dims->N = opts->N2;
}



int ocp_qp_partial_condensing_memory_calculate_size(ocp_qp_dims *dims, ocp_qp_partial_condensing_opts *opts)
{
    int size = 0;

    // populate dimensions of new ocp_qp based on N2
    opts->pcond_dims->N = opts->N2;
    d_compute_qp_dim_ocp2ocp(dims, opts->pcond_dims);

    size += sizeof(ocp_qp_partial_condensing_memory);
    size += sizeof(struct d_cond_qp_ocp2ocp_workspace);
    size += d_memsize_cond_qp_ocp2ocp(dims, opts->pcond_dims);

    return size;
}



ocp_qp_partial_condensing_memory *ocp_qp_partial_condensing_memory_assign(ocp_qp_dims *dims,
    ocp_qp_partial_condensing_opts *opts, void *raw_memory)
{
    char *c_ptr = (char *)raw_memory;

    ocp_qp_partial_condensing_memory *mem = (ocp_qp_partial_condensing_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_partial_condensing_memory);
    //
    mem->hpipm_workspace = (struct d_cond_qp_ocp2ocp_workspace *)c_ptr;
    c_ptr += sizeof(struct d_cond_qp_ocp2ocp_workspace);

    // hpipm workspace structure
    assert((size_t)c_ptr % 8 == 0 && "double not 8-byte aligned!");

    d_create_cond_qp_ocp2ocp(dims, opts->pcond_dims, mem->hpipm_workspace, c_ptr);
    c_ptr += mem->hpipm_workspace->memsize;

    mem->qp_in = NULL;  // initialized when partial condensing routine is called

    assert((char*)raw_memory + ocp_qp_partial_condensing_memory_calculate_size(dims, opts) == c_ptr);

    return mem;
}



int ocp_qp_partial_condensing_workspace_calculate_size(ocp_qp_dims *dims, ocp_qp_partial_condensing_opts *opts)
{
    return 0;
}



void ocp_qp_partial_condensing(ocp_qp_in *in, ocp_qp_in *out, ocp_qp_partial_condensing_opts *opts,
    ocp_qp_partial_condensing_memory *mem, void *work)
{
    // save pointers to ocp_qp_in in memory (needed for expansion)
    mem->qp_in = in;
    mem->pcond_qp_in = out;

    // convert to partially condensed qp structure
    d_cond_qp_ocp2ocp(in, out, mem->hpipm_workspace);
}



void ocp_qp_partial_expansion(ocp_qp_out *in, ocp_qp_out *out, ocp_qp_partial_condensing_opts *opts,
    ocp_qp_partial_condensing_memory *mem, void *work)
{
    d_expand_sol_ocp2ocp(mem->qp_in, mem->pcond_qp_in, in, out, mem->hpipm_workspace);
}
