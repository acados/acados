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

// hpipm
#include "hpipm_d_dense_qp.h"
#include "hpipm_d_dense_qp_sol.h"
#include "hpipm_d_dense_qp_ipm.h"
// acados
#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/mem.h"


int dense_qp_hpipm_calculate_args_size(dense_qp_in *qp_in) {

    int size = 0;
    size += sizeof(dense_qp_hpipm_args);
    size += sizeof(struct d_dense_qp_ipm_arg);
    size += d_memsize_dense_qp_ipm_arg(qp_in);

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



char *dense_qp_hpipm_assign_args(dense_qp_in *qp_in, dense_qp_hpipm_args **args, void *mem) {

    char *c_ptr = (char *) mem;

    *args = (dense_qp_hpipm_args *) c_ptr;
    c_ptr += sizeof(dense_qp_hpipm_args);

    (*args)->hpipm_args = (struct d_dense_qp_ipm_arg *) c_ptr;
    c_ptr += sizeof(struct d_dense_qp_ipm_arg);

    align_char_to(8, &c_ptr);
    d_create_dense_qp_ipm_arg(qp_in, (*args)->hpipm_args, c_ptr);
    c_ptr += d_memsize_dense_qp_ipm_arg(qp_in);

    return c_ptr;
}



void dense_qp_hpipm_initialize_default_args(dense_qp_hpipm_args *args) {

    d_set_default_dense_qp_ipm_arg(args->hpipm_args);
    // overwrite some default options
    args->hpipm_args->res_g_max = 1e-6;
    args->hpipm_args->res_b_max = 1e-8;
    args->hpipm_args->res_d_max = 1e-8;
    args->hpipm_args->res_m_max = 1e-8;
    args->hpipm_args->iter_max = 50;
    args->hpipm_args->stat_max = 50;
    args->hpipm_args->alpha_min = 1e-8;
    args->hpipm_args->mu0 = 1;
}



int dense_qp_hpipm_calculate_memory_size(dense_qp_in *qp_in, dense_qp_hpipm_args *args) {

    int size = 0;
    size += sizeof(dense_qp_hpipm_memory);

    size += sizeof(struct d_dense_qp_ipm_workspace);  // ipm_workspace

    size += d_memsize_dense_qp_ipm(qp_in, args->hpipm_args);

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



char *dense_qp_hpipm_assign_memory(dense_qp_in *qp_in, dense_qp_hpipm_args *args,
        void **mem_, void *raw_memory) {

    dense_qp_hpipm_memory **hpipm_memory = (dense_qp_hpipm_memory **) mem_;

    // char pointer
    char *c_ptr = (char *)raw_memory;

    *hpipm_memory = (dense_qp_hpipm_memory *) c_ptr;
    c_ptr += sizeof(dense_qp_hpipm_memory);

    //
    (*hpipm_memory)->hpipm_workspace = (struct d_dense_qp_ipm_workspace *)c_ptr;
    c_ptr += sizeof(struct d_dense_qp_ipm_workspace);

    struct d_dense_qp_ipm_workspace *ipm_workspace = (*hpipm_memory)->hpipm_workspace;

    // ipm workspace structure
    align_char_to(8, &c_ptr);
    d_create_dense_qp_ipm(qp_in, args->hpipm_args, ipm_workspace, c_ptr);
    c_ptr += ipm_workspace->memsize;

    return c_ptr;
}



int dense_qp_hpipm(dense_qp_in *qp_in, dense_qp_out *qp_out, void *args_, void *mem_) {

    dense_qp_hpipm_args *args = (dense_qp_hpipm_args *) args_;
    dense_qp_hpipm_memory *memory = (dense_qp_hpipm_memory *) mem_;

    // initialize return code
    int acados_status = ACADOS_SUCCESS;
    // solve ipm
    int hpipm_status = d_solve_dense_qp_ipm(qp_in, qp_out, args->hpipm_args, memory->hpipm_workspace);

    // check max number of iterations
    // TODO(dimitris): check ACADOS_MIN_STEP (not implemented in HPIPM yet)
    if (hpipm_status == 1) acados_status = ACADOS_MAXITER;

    // return
    return acados_status;
}
