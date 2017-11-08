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
#if defined(RUNTIME_CHECKS)
#include <assert.h>
#endif
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



void dummy_dense_qp_in(dense_qp_in *qpd_in, ocp_qp_dims *dims) {

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



int ocp_qp_condensing_calculate_args_size(ocp_qp_dims *dims) {

    int size = 0;
    size += sizeof(ocp_qp_condensing_args);
    return size;
}



ocp_qp_condensing_args *ocp_qp_condensing_assign_args(ocp_qp_dims *dims, void *mem) {

    ocp_qp_condensing_args *args;

    char *c_ptr = (char *) mem;

    args = (ocp_qp_condensing_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_args);

#if defined(RUNTIME_CHECKS)
    assert((char*)mem + ocp_qp_condensing_calculate_args_size(dims) >= c_ptr);
#endif
    return args;
}



int ocp_qp_condensing_calculate_memory_size(ocp_qp_dims *dims, ocp_qp_condensing_args *args) {

    int size = sizeof(ocp_qp_condensing_memory);

    size += sizeof(struct d_cond_qp_ocp2dense_workspace);

    // // dummy dense qp
    // dense_qp_in qpd_in;
    // dummy_dense_qp_in(&qpd_in, qp_in);

    size += d_memsize_cond_qp_ocp2dense(dims);

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}


char *assign_ocp_qp_condensing_memory(ocp_qp_in *in, ocp_qp_condensing_memory **memory,
    void *raw_memory) {

    // char pointer
    char *c_ptr = (char *)raw_memory;

    *memory = (ocp_qp_condensing_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_memory);
    // pointer to ocp_qp, needed for the expansion (not copied to memory)
    (*memory)->qp_in = in;
    //
    (*memory)->hpipm_workspace = (struct d_cond_qp_ocp2dense_workspace *)c_ptr;
    c_ptr += sizeof(struct d_cond_qp_ocp2dense_workspace);

    // hpipm workspace structure
    align_char_to(8, &c_ptr);
    struct d_cond_qp_ocp2dense_workspace *hpipm_workspace = (*memory)->hpipm_workspace;
    d_create_cond_qp_ocp2dense(in->size, hpipm_workspace, c_ptr);
    c_ptr += hpipm_workspace->memsize;

    return c_ptr;
}



void ocp_qp_condensing(ocp_qp_in *in, dense_qp_in *out, ocp_qp_condensing_args *args,
    ocp_qp_condensing_memory *mem) {

    // convert to dense qp structure
    d_cond_qp_ocp2dense(in, out, mem->hpipm_workspace);
}



// TODO(dimitris): Remove qp from inputs
void ocp_qp_expansion(dense_qp_out *in, ocp_qp_out *out, ocp_qp_condensing_args *args,
    ocp_qp_condensing_memory *mem) {

    d_expand_sol_dense2ocp(mem->qp_in, in, out, mem->hpipm_workspace);
}
