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
#include <stdio.h>
// acados
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/ocp_qp/ocp_qp_condensing.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/dense_qp/dense_qp_qpoases.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"


int ocp_qp_condensing_qpoases_calculate_args_size(ocp_qp_dims *dims) {

    int size = 0;
    size += sizeof(ocp_qp_condensing_qpoases_args);

    // dummy dense qp
    dense_qp_in qpd_in;
    dummy_dense_qp_in(&qpd_in, dims);

    size += sizeof(ocp_qp_condensing_args);
    size += sizeof(dense_qp_qpoases_args);

    size += dense_qp_qpoases_calculate_args_size(&qpd_in);
    size += ocp_qp_condensing_calculate_args_size(dims);

    return size;
}



void *ocp_qp_condensing_qpoases_assign_args(ocp_qp_dims *dims, void *mem)
{
    ocp_qp_condensing_qpoases_args *args;
    char *c_ptr = (char *) mem;

    args = (ocp_qp_condensing_qpoases_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_qpoases_args);

    args->cond_args = (ocp_qp_condensing_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_args);

    args->solver_args = (dense_qp_qpoases_args *) c_ptr;
    c_ptr += sizeof(dense_qp_qpoases_args);

    // dummy dense qp
    dense_qp_in qpd_in;
    dummy_dense_qp_in(&qpd_in, dims);

    args->solver_args = dense_qp_qpoases_assign_args(&qpd_in, c_ptr);
    c_ptr += dense_qp_qpoases_calculate_args_size(&qpd_in);

    args->cond_args = ocp_qp_condensing_assign_args(dims, c_ptr);
    c_ptr += ocp_qp_condensing_calculate_args_size(dims);

#if defined(RUNTIME_CHECKS)
    assert((char*)mem + ocp_qp_condensing_qpoases_calculate_args_size(dims) >= c_ptr);
#endif
    return (void*)args;
}



void ocp_qp_condensing_qpoases_initialize_default_args(void *args_)
{
    ocp_qp_condensing_qpoases_args *args = (ocp_qp_condensing_qpoases_args *)args_;
    args->solver_args->max_cputime = 1000.0;
    args->solver_args->warm_start = 0;
    args->solver_args->max_nwsr = 1000;
}



int ocp_qp_condensing_qpoases_calculate_memory_size(ocp_qp_dims *dims,
    ocp_qp_condensing_qpoases_args *args) {

    int size = 0;
    size += sizeof(ocp_qp_condensing_qpoases_memory);

    size += sizeof(ocp_qp_condensing_memory);
    size += sizeof(dense_qp_qpoases_memory);
    size += sizeof(dense_qp_in);
    size += sizeof(dense_qp_out);

    // dummy dense qp
    dense_qp_in qpd_in;
    dummy_dense_qp_in(&qpd_in, dims);

    size += ocp_qp_condensing_calculate_memory_size(dims, args->cond_args);
    size += dense_qp_qpoases_calculate_memory_size(&qpd_in, args->solver_args);
    size += dense_qp_in_calculate_size(qpd_in.nv, qpd_in.ne, qpd_in.nb, qpd_in.ng, qpd_in.ns);
    size += dense_qp_out_calculate_size(qpd_in.nv, qpd_in.ne, qpd_in.nb, qpd_in.ng, qpd_in.ns);

    make_int_multiple_of(8, &size);
    size += 4 * 8;

    return size;
}



char *ocp_qp_condensing_qpoases_assign_memory(ocp_qp_dims *dims, ocp_qp_condensing_qpoases_args *args,
        void **mem_, void *raw_memory) {

    ocp_qp_condensing_qpoases_memory **qpoases_memory = (ocp_qp_condensing_qpoases_memory **) mem_;

    // char pointer
    char *c_ptr = (char *)raw_memory;

    *qpoases_memory = (ocp_qp_condensing_qpoases_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_qpoases_memory);

    //
    (*qpoases_memory)->condensing_memory = (ocp_qp_condensing_memory *)c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_memory);
    //
    (*qpoases_memory)->solver_memory = (dense_qp_qpoases_memory *)c_ptr;
    c_ptr += sizeof(dense_qp_qpoases_memory);
    //
    (*qpoases_memory)->qpd_in = (dense_qp_in *)c_ptr;
    c_ptr += sizeof(dense_qp_in);
    //
    (*qpoases_memory)->qpd_out = (dense_qp_out *)c_ptr;
    c_ptr += sizeof(dense_qp_out);

    // dummy dense qp
    dense_qp_in qpd_in;
    dummy_dense_qp_in(&qpd_in, dims);

    align_char_to(8, &c_ptr);
    assign_ocp_qp_condensing_memory(dims, &(*qpoases_memory)->condensing_memory, c_ptr);
    c_ptr += ocp_qp_condensing_calculate_memory_size(dims, args->cond_args);

    // TODO(dimitris): TRANSFER ASSERT INSIDE ASSIGN RATHER THAN CREATE! ASSERTS OF DENSE QP NOT CHECKED
    align_char_to(8, &c_ptr);
    dense_qp_qpoases_assign_memory(&qpd_in, args->solver_args, (void**)&(*qpoases_memory)->solver_memory, c_ptr);
    c_ptr += dense_qp_qpoases_calculate_memory_size(&qpd_in, args->solver_args);

    align_char_to(8, &c_ptr);
    assign_dense_qp_in(qpd_in.nv, qpd_in.ne, qpd_in.nb, qpd_in.ng, qpd_in.ns, &(*qpoases_memory)->qpd_in, c_ptr);
    c_ptr += dense_qp_in_calculate_size(qpd_in.nv, qpd_in.ne, qpd_in.nb, qpd_in.ng, qpd_in.ns);

    align_char_to(8, &c_ptr);
    assign_dense_qp_out(qpd_in.nv, qpd_in.ne, qpd_in.nb, qpd_in.ng, qpd_in.ns, &(*qpoases_memory)->qpd_out, c_ptr);
    c_ptr += dense_qp_out_calculate_size(qpd_in.nv, qpd_in.ne, qpd_in.nb, qpd_in.ng, qpd_in.ns);

    return c_ptr;
}



int ocp_qp_condensing_qpoases(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *mem_) {

    ocp_qp_condensing_qpoases_args *args = (ocp_qp_condensing_qpoases_args *) args_;
    ocp_qp_condensing_qpoases_memory *memory = (ocp_qp_condensing_qpoases_memory *) mem_;

    // initialize return code
    int acados_status = ACADOS_SUCCESS;

    // condense
    ocp_qp_condensing(qp_in, memory->qpd_in, args->cond_args, memory->condensing_memory);

    // solve ipm
    int qpoases_status = dense_qp_qpoases(memory->qpd_in, memory->qpd_out, args->solver_args, memory->solver_memory);

    // expand
    ocp_qp_expansion(memory->qpd_out, qp_out, args->cond_args, memory->condensing_memory);

    // check max number of iterations
    // TODO(dimitris): check ACADOS_MIN_STEP (not implemented in HPIPM yet)
    if (qpoases_status == 1) acados_status = ACADOS_MAXITER;

    // return
    return acados_status;
}
