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
// hpipm
// #include "hpipm_d_ocp_qp.h"
// #include "hpipm_d_ocp_qp_sol.h"
// #include "hpipm_d_dense_qp.h"
// #include "hpipm_d_dense_qp_sol.h"
// #include "hpipm_d_cond.h"
// #include "hpipm_d_dense_qp_ipm.h"
// acados
#include "acados/ocp_qp/ocp_qp_condensing_hpipm.h"
#include "acados/ocp_qp/ocp_qp_condensing.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"


int ocp_qp_condensing_hpipm_calculate_args_size(ocp_qp_dims *dims) {

    int size = 0;
    size += sizeof(ocp_qp_condensing_hpipm_args);

    // dummy dense qp
    dense_qp_in qpd_in;
    dummy_dense_qp_in(&qpd_in, dims);

    size += dense_qp_hpipm_calculate_args_size(&qpd_in);
    size += ocp_qp_condensing_calculate_args_size(dims);

    return size;
}



void *ocp_qp_condensing_hpipm_assign_args(ocp_qp_dims *dims, void *mem)
{
    char *c_ptr = (char *) mem;

    ocp_qp_condensing_hpipm_args *args = (ocp_qp_condensing_hpipm_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_hpipm_args);

    // dummy dense qp
    dense_qp_in qpd_in;
    dummy_dense_qp_in(&qpd_in, dims);

    args->solver_args = dense_qp_hpipm_assign_args(&qpd_in, c_ptr);
    c_ptr += dense_qp_hpipm_calculate_args_size(&qpd_in);

    args->cond_args = ocp_qp_condensing_assign_args(dims, c_ptr);
    c_ptr += ocp_qp_condensing_calculate_args_size(dims);

#if defined(RUNTIME_CHECKS)
    assert((char*)mem + ocp_qp_condensing_hpipm_calculate_args_size(dims) >= c_ptr);
#endif
    return (void *)args;
}



void ocp_qp_condensing_hpipm_initialize_default_args(void *args_)
{
    ocp_qp_condensing_hpipm_args *args = (ocp_qp_condensing_hpipm_args *)args_;

    d_set_default_dense_qp_ipm_arg(args->solver_args->hpipm_args);
    // overwrite some default options
    args->solver_args->hpipm_args->res_g_max = 1e-6;
    args->solver_args->hpipm_args->res_b_max = 1e-8;
    args->solver_args->hpipm_args->res_d_max = 1e-8;
    args->solver_args->hpipm_args->res_m_max = 1e-8;
    args->solver_args->hpipm_args->iter_max = 50;
    args->solver_args->hpipm_args->stat_max = 50;
    args->solver_args->hpipm_args->alpha_min = 1e-8;
    args->solver_args->hpipm_args->mu0 = 1;
}



int ocp_qp_condensing_hpipm_calculate_memory_size(ocp_qp_dims *dims, void *args_)
{
    ocp_qp_condensing_hpipm_args *args = (ocp_qp_condensing_hpipm_args *)args_;

    int size = 0;
    size += sizeof(ocp_qp_condensing_hpipm_memory);

    // dummy dense qp
    dense_qp_in qpd_in;
    dummy_dense_qp_in(&qpd_in, dims);

    size += ocp_qp_condensing_calculate_memory_size(dims, args->cond_args);
    size += dense_qp_hpipm_calculate_memory_size(&qpd_in, args->solver_args);
    size += dense_qp_in_calculate_size(qpd_in.nv, qpd_in.ne, qpd_in.nb, qpd_in.ng, qpd_in.ns);
    size += dense_qp_out_calculate_size(qpd_in.nv, qpd_in.ne, qpd_in.nb, qpd_in.ng, qpd_in.ns);

    make_int_multiple_of(8, &size);
    size += 4 * 8;

    return size;
}



void *ocp_qp_condensing_hpipm_assign_memory(ocp_qp_dims *dims, void *args_, void *raw_memory)
{
    ocp_qp_condensing_hpipm_args *args = (ocp_qp_condensing_hpipm_args *)args_;

    char *c_ptr = (char *)raw_memory;

    ocp_qp_condensing_hpipm_memory *mem = (ocp_qp_condensing_hpipm_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_hpipm_memory);

    // dummy dense qp
    dense_qp_in qpd_in;
    dummy_dense_qp_in(&qpd_in, dims);

    align_char_to(8, &c_ptr);
    mem->condensing_memory = assign_ocp_qp_condensing_memory(dims, args->cond_args, c_ptr);
    c_ptr += ocp_qp_condensing_calculate_memory_size(dims, args->cond_args);

    align_char_to(8, &c_ptr);
    mem->solver_memory = dense_qp_hpipm_assign_memory(&qpd_in, args->solver_args, c_ptr);
    c_ptr += dense_qp_hpipm_calculate_memory_size(&qpd_in, args->solver_args);

    align_char_to(8, &c_ptr);
    assign_dense_qp_in(qpd_in.nv, qpd_in.ne, qpd_in.nb, qpd_in.ng, qpd_in.ns, &mem->qpd_in, c_ptr);
    c_ptr += dense_qp_in_calculate_size(qpd_in.nv, qpd_in.ne, qpd_in.nb, qpd_in.ng, qpd_in.ns);

    align_char_to(8, &c_ptr);
    assign_dense_qp_out(qpd_in.nv, qpd_in.ne, qpd_in.nb, qpd_in.ng, qpd_in.ns, &mem->qpd_out, c_ptr);
    c_ptr += dense_qp_out_calculate_size(qpd_in.nv, qpd_in.ne, qpd_in.nb, qpd_in.ng, qpd_in.ns);

#if defined(RUNTIME_CHECKS)
    assert((char *) raw_memory + ocp_qp_condensing_hpipm_calculate_memory_size(dims, args_) >= c_ptr);
#endif
    return mem;
}



int ocp_qp_condensing_hpipm(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *mem_) {

    ocp_qp_condensing_hpipm_args *args = (ocp_qp_condensing_hpipm_args *) args_;
    ocp_qp_condensing_hpipm_memory *memory = (ocp_qp_condensing_hpipm_memory *) mem_;

    // initialize return code
    int acados_status = ACADOS_SUCCESS;

    // condense
    ocp_qp_condensing(qp_in, memory->qpd_in, args->cond_args, memory->condensing_memory);

    // solve ipm
    int hpipm_status = dense_qp_hpipm(memory->qpd_in, memory->qpd_out, args->solver_args, memory->solver_memory);

    // expand
    ocp_qp_expansion(memory->qpd_out, qp_out, args->cond_args, memory->condensing_memory);

    // check max number of iterations
    // TODO(dimitris): check ACADOS_MIN_STEP (not implemented in HPIPM yet)
    if (hpipm_status == 1) acados_status = ACADOS_MAXITER;

    // return
    return acados_status;
}
