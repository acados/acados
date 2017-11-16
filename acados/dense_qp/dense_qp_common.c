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
#include "hpipm_d_dense_qp.h"
#include "hpipm_d_dense_qp_sol.h"
// acados
#include "acados/utils/types.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/dense_qp/dense_qp_qpoases.h"



int dense_qp_in_calculate_size(dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_in);
    size += sizeof(dense_qp_dims);
    size += d_memsize_dense_qp(dims);
    return size;
}



dense_qp_in *assign_dense_qp_in(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *)raw_memory;

    dense_qp_in *qp_in = (dense_qp_in *) c_ptr;
    c_ptr += sizeof(dense_qp_in);

    d_create_dense_qp(dims, qp_in, c_ptr);
    c_ptr += d_memsize_dense_qp(dims);

    qp_in->dim = (dense_qp_dims *) c_ptr;
    c_ptr += sizeof(dense_qp_dims);

    qp_in->dim->nv = dims->nv;
    qp_in->dim->ne = dims->ne;
    qp_in->dim->nb = dims->nb;
    qp_in->dim->ng = dims->ng;
    qp_in->dim->ns = dims->ns;

    assert((char*)raw_memory + dense_qp_in_calculate_size(dims) == c_ptr);

    return qp_in;
}



int dense_qp_out_calculate_size(dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_out);
    size += d_memsize_dense_qp_sol(dims);

    return size;
}



dense_qp_out *assign_dense_qp_out(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_out *qp_out = (dense_qp_out *) c_ptr;
    c_ptr += sizeof(dense_qp_out);

    d_create_dense_qp_sol(dims, qp_out, c_ptr);
    c_ptr += d_memsize_dense_qp_sol(dims);

    assert((char*)raw_memory + dense_qp_out_calculate_size(dims) >= c_ptr);

    return qp_out;
}



void set_dense_qp_solver_fun_ptrs(dense_qp_solver_t qp_solver_name, dense_qp_solver *qp_solver)
{
    switch (qp_solver_name)
    {
        case DENSE_QP_HPIPM:
            qp_solver->calculate_args_size = &dense_qp_hpipm_calculate_args_size;
            qp_solver->assign_args = &dense_qp_hpipm_assign_args;
            qp_solver->initialize_default_args = &dense_qp_hpipm_initialize_default_args;
            qp_solver->calculate_memory_size = &dense_qp_hpipm_calculate_memory_size;
            qp_solver->assign_memory = &dense_qp_hpipm_assign_memory;
            qp_solver->fun = &dense_qp_hpipm;
            break;
        case DENSE_QP_QPOASES:
            qp_solver->calculate_args_size = &dense_qp_qpoases_calculate_args_size;
            qp_solver->assign_args = &dense_qp_qpoases_assign_args;
            qp_solver->initialize_default_args = &dense_qp_qpoases_initialize_default_args;
            qp_solver->calculate_memory_size = &dense_qp_qpoases_calculate_memory_size;
            qp_solver->assign_memory = &dense_qp_qpoases_assign_memory;
            qp_solver->fun = &dense_qp_qpoases;
            break;
        default:
            qp_solver->calculate_args_size = NULL;
            qp_solver->assign_args = NULL;
            qp_solver->initialize_default_args = NULL;
            qp_solver->calculate_memory_size = NULL;
            qp_solver->assign_memory = NULL;
            qp_solver->fun = NULL;
    }
}
