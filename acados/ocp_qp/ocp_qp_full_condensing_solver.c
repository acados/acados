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
#include <stdio.h>
// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/ocp_qp/ocp_qp_full_condensing_solver.h"
#include "acados/ocp_qp/ocp_qp_full_condensing.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"


int ocp_qp_full_condensing_solver_calculate_args_size(ocp_qp_dims *dims, void *solver_)
{
    dense_qp_solver *solver = (dense_qp_solver *)solver_;

    int size = 0;
    size += sizeof(ocp_qp_full_condensing_solver_args);
    size += sizeof(dense_qp_solver);

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    size += ocp_qp_full_condensing_calculate_args_size(dims);
    size += solver->calculate_args_size(&ddims);

    return size;
}



void *ocp_qp_full_condensing_solver_assign_args(ocp_qp_dims *dims, void *solver_, void *raw_memory)
{
    dense_qp_solver *solver = (dense_qp_solver *)solver_;

    char *c_ptr = (char *) raw_memory;

    ocp_qp_full_condensing_solver_args *args = (ocp_qp_full_condensing_solver_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_full_condensing_solver_args);

    args->solver = (dense_qp_solver*) c_ptr;
    c_ptr += sizeof(dense_qp_solver);

    copy_module_pointers_to_args(args->solver, solver);

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    args->cond_args = ocp_qp_full_condensing_assign_args(dims, c_ptr);
    c_ptr += ocp_qp_full_condensing_calculate_args_size(dims);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    args->solver_args = args->solver->assign_args(&ddims, c_ptr);
    c_ptr += args->solver->calculate_args_size(&ddims);

    assert((char*)raw_memory + ocp_qp_full_condensing_solver_calculate_args_size(dims, solver) == c_ptr);

    return (void*)args;
}



void ocp_qp_full_condensing_solver_initialize_default_args(void *args_)
{
    ocp_qp_full_condensing_solver_args *args = (ocp_qp_full_condensing_solver_args *)args_;
    // ocp_qp_full_condensing_initialize_default_args(args->cond_args);
    args->solver->initialize_default_args(args->solver_args);
}



int ocp_qp_full_condensing_solver_calculate_memory_size(ocp_qp_dims *dims, void *args_)
{
    ocp_qp_full_condensing_solver_args *args = (ocp_qp_full_condensing_solver_args *)args_;

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    int size = 0;
    size += sizeof(ocp_qp_full_condensing_solver_memory);

    size += ocp_qp_full_condensing_calculate_memory_size(dims, args->cond_args);
    size += args->solver->calculate_memory_size(&ddims, args->solver_args);

    size += dense_qp_in_calculate_size(&ddims);
    size += dense_qp_out_calculate_size(&ddims);

    return size;
}



void *ocp_qp_full_condensing_solver_assign_memory(ocp_qp_dims *dims, void *args_, void *raw_memory)
{
    ocp_qp_full_condensing_solver_args *args = (ocp_qp_full_condensing_solver_args *)args_;

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    char *c_ptr = (char *)raw_memory;

    ocp_qp_full_condensing_solver_memory *mem = (ocp_qp_full_condensing_solver_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_full_condensing_solver_memory);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    mem->cond_memory = ocp_qp_full_condensing_assign_memory(dims, args->cond_args, c_ptr);
    c_ptr += ocp_qp_full_condensing_calculate_memory_size(dims, args->cond_args);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    mem->solver_memory = args->solver->assign_memory(&ddims, args->solver_args, c_ptr);
    c_ptr += args->solver->calculate_memory_size(&ddims, args->solver_args);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    mem->qpd_in = assign_dense_qp_in(&ddims, c_ptr);
    c_ptr += dense_qp_in_calculate_size(&ddims);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    mem->qpd_out = assign_dense_qp_out(&ddims, c_ptr);
    c_ptr += dense_qp_out_calculate_size(&ddims);

    assert((char *) raw_memory + ocp_qp_full_condensing_solver_calculate_memory_size(dims, args_) == c_ptr);

    return mem;
}



int ocp_qp_full_condensing_solver_calculate_workspace_size(ocp_qp_dims *dims, void *args_)
{
    return sizeof(ocp_qp_full_condensing_solver_workspace);
}



int ocp_qp_full_condensing_solver(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *mem_, void *work_)
{
    ocp_qp_full_condensing_solver_args *args = (ocp_qp_full_condensing_solver_args *) args_;
    ocp_qp_full_condensing_solver_memory *memory = (ocp_qp_full_condensing_solver_memory *) mem_;
    // TODO(dimitris): also assign workspace once it contains data that need to be casted..
    ocp_qp_full_condensing_solver_workspace *work = (ocp_qp_full_condensing_solver_workspace *) work_;

    // condense
    ocp_qp_full_condensing(qp_in, memory->qpd_in, args->cond_args, memory->cond_memory, work->cond_work);

    // solve qp
    int solver_status = args->solver->fun(memory->qpd_in, memory->qpd_out, args->solver_args, memory->solver_memory, work->solver_workspace);

    // expand
    ocp_qp_expansion(memory->qpd_out, qp_out, args->cond_args, memory->cond_memory, work->cond_work);

    // return
    return solver_status;
}
