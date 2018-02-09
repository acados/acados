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
// acados
#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"



int dense_qp_hpipm_calculate_args_size(dense_qp_dims *dims, void *submodules_)
{
    int size = 0;
    size += sizeof(dense_qp_hpipm_args);
    size += sizeof(struct d_dense_qp_ipm_arg);
    size += d_memsize_dense_qp_ipm_arg(dims);

    return size;
}



void *dense_qp_hpipm_assign_args(dense_qp_dims *dims, void **submodules_, void *raw_memory)
{
    dense_qp_hpipm_args *args;

    char *c_ptr = (char *) raw_memory;

    args = (dense_qp_hpipm_args *) c_ptr;
    c_ptr += sizeof(dense_qp_hpipm_args);

    args->hpipm_args = (struct d_dense_qp_ipm_arg *) c_ptr;
    c_ptr += sizeof(struct d_dense_qp_ipm_arg);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    d_create_dense_qp_ipm_arg(dims, args->hpipm_args, c_ptr);
    c_ptr += d_memsize_dense_qp_ipm_arg(dims);

    assert((char*)raw_memory + dense_qp_hpipm_calculate_args_size(dims, *submodules_) == c_ptr);

    // Update submodules pointer
    *submodules_ = NULL;

    return (void *)args;
}



void *dense_qp_hpipm_copy_args(dense_qp_dims *dims, void *raw_memory, void *source_)
{
    dense_qp_hpipm_args *source = (dense_qp_hpipm_args *)source_;
    dense_qp_hpipm_args *dest;

    void *submodules;

    dest = (dense_qp_hpipm_args *) dense_qp_hpipm_assign_args(dims, &submodules, raw_memory);

    *dest->hpipm_args = *source->hpipm_args;

    return (void *)dest;
}



void dense_qp_hpipm_initialize_default_args(void *args_)
{
    dense_qp_hpipm_args *args = (dense_qp_hpipm_args *)args_;

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



int dense_qp_hpipm_calculate_memory_size(dense_qp_dims *dims, void *args_)
{
    dense_qp_hpipm_args *args = (dense_qp_hpipm_args *)args_;

    int size = 0;
    size += sizeof(dense_qp_hpipm_memory);
    size += sizeof(struct d_dense_qp_ipm_workspace);

    size += d_memsize_dense_qp_ipm(dims, args->hpipm_args);

    return size;
}



void *dense_qp_hpipm_assign_memory(dense_qp_dims *dims, void *args_, void *raw_memory)
{
    dense_qp_hpipm_args *args = (dense_qp_hpipm_args *)args_;
    dense_qp_hpipm_memory *mem;

    char *c_ptr = (char *)raw_memory;

    mem = (dense_qp_hpipm_memory *) c_ptr;
    c_ptr += sizeof(dense_qp_hpipm_memory);

    mem->hpipm_workspace = (struct d_dense_qp_ipm_workspace *)c_ptr;
    c_ptr += sizeof(struct d_dense_qp_ipm_workspace);

    struct d_dense_qp_ipm_workspace *ipm_workspace = mem->hpipm_workspace;

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    // ipm workspace structure
    d_create_dense_qp_ipm(dims, args->hpipm_args, ipm_workspace, c_ptr);
    c_ptr += ipm_workspace->memsize;

    assert((char *)raw_memory + dense_qp_hpipm_calculate_memory_size(dims, args) == c_ptr);

    return mem;
}



int dense_qp_hpipm_calculate_workspace_size(dense_qp_dims *dims, void *args_)
{
    return 0;
}



int dense_qp_hpipm(dense_qp_in *qp_in, dense_qp_out *qp_out, void *args_, void *mem_, void *work_)
{
    dense_qp_info *info = (dense_qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer;

    acados_tic(&tot_timer);

    // cast structures
    dense_qp_hpipm_args *args = (dense_qp_hpipm_args *) args_;
    dense_qp_hpipm_memory *memory = (dense_qp_hpipm_memory *) mem_;

    // solve ipm
    acados_tic(&qp_timer);
    int hpipm_status = d_solve_dense_qp_ipm(qp_in, qp_out, args->hpipm_args, memory->hpipm_workspace);

    info->solve_QP_time = acados_toc(&qp_timer);
    info->interface_time = 0;  // there are no conversions for hpipm
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = memory->hpipm_workspace->iter;

    // check exit conditions
    int acados_status = hpipm_status;
    if (hpipm_status == 0) acados_status = ACADOS_SUCCESS;
    if (hpipm_status == 1) acados_status = ACADOS_MAXITER;
    if (hpipm_status == 2) acados_status = ACADOS_MINSTEP;
    return acados_status;
}
