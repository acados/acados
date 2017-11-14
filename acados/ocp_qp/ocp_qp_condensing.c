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
// hpipm
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_dense_qp_ipm.h"
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_cond.h"
// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/ocp_qp/ocp_qp_condensing.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"


// TODO(dimitris): implement dense_dims in hpipm and get rid of this file
void dummy_dense_qp_in(dense_qp_in *qpd_in, ocp_qp_dims *dims)
{
    // compute dense qp size
    int nvd = 0;
    int ned = 0;
    int nbd = 0;
    int ngd = 0;
    int nsd = 0;

    d_compute_qp_size_ocp2dense(dims, &nvd, &ned, &nbd, &ngd, &nsd);

    // dummy dense qp
    qpd_in->nv = nvd;
    qpd_in->ne = ned;
    qpd_in->nb = nbd;
    qpd_in->ng = ngd;
    qpd_in->ns = nsd;
}



int ocp_qp_condensing_calculate_args_size(ocp_qp_dims *dims)
{
    int size = 0;
    size += sizeof(ocp_qp_condensing_args);
    return size;
}



ocp_qp_condensing_args *ocp_qp_condensing_assign_args(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_condensing_args *args = (ocp_qp_condensing_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_args);

    assert((char*)raw_memory + ocp_qp_condensing_calculate_args_size(dims) >= c_ptr);

    return args;
}



int ocp_qp_condensing_calculate_memory_size(ocp_qp_dims *dims, ocp_qp_condensing_args *args)
{
    int size = 0;

    size += sizeof(ocp_qp_condensing_memory);
    size += sizeof(struct d_cond_qp_ocp2dense_workspace);
    size += d_memsize_cond_qp_ocp2dense(dims);

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



void *assign_ocp_qp_condensing_memory(ocp_qp_dims *dims, ocp_qp_condensing_args *args, void *raw_memory)
{
    char *c_ptr = (char *)raw_memory;

    ocp_qp_condensing_memory *mem = (ocp_qp_condensing_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_memory);
    //
    mem->hpipm_workspace = (struct d_cond_qp_ocp2dense_workspace *)c_ptr;
    c_ptr += sizeof(struct d_cond_qp_ocp2dense_workspace);

    // hpipm workspace structure
    align_char_to(8, &c_ptr);
    d_create_cond_qp_ocp2dense(dims, mem->hpipm_workspace, c_ptr);
    c_ptr += mem->hpipm_workspace->memsize;

    assert((char*)raw_memory + ocp_qp_condensing_calculate_memory_size(dims, args) >= c_ptr);

    return mem;
}



void ocp_qp_condensing(ocp_qp_in *in, dense_qp_in *out, ocp_qp_condensing_args *args, ocp_qp_condensing_memory *mem)
{
    // save pointer to ocp_qp_in to memory (needed for expansion)
    mem->qp_in = in;

    // convert to dense qp structure
    d_cond_qp_ocp2dense(in, out, mem->hpipm_workspace);
}



void ocp_qp_expansion(dense_qp_out *in, ocp_qp_out *out, ocp_qp_condensing_args *args, ocp_qp_condensing_memory *mem)
{
    d_expand_sol_dense2ocp(mem->qp_in, in, out, mem->hpipm_workspace);
}
