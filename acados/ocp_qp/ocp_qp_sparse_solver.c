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
#include "acados/ocp_qp/ocp_qp_sparse_solver.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"
#include "acados/utils/timing.h"
#include "acados/utils/mem.h"


int ocp_qp_sparse_solver_calculate_args_size(ocp_qp_dims *dims, void *solver_)
{
    ocp_qp_solver *solver = (ocp_qp_solver *)solver_;

    int size = 0;
    size += sizeof(ocp_qp_sparse_solver_args);
    size += sizeof(ocp_qp_solver);

    size += ocp_qp_partial_condensing_calculate_args_size(dims);
    size += solver->calculate_args_size(dims);

    return size;
}



void *ocp_qp_sparse_solver_assign_args(ocp_qp_dims *dims, void *solver_, void *raw_memory)
{
    ocp_qp_solver *solver = (ocp_qp_solver *)solver_;

    char *c_ptr = (char *) raw_memory;

    ocp_qp_sparse_solver_args *args = (ocp_qp_sparse_solver_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_sparse_solver_args);

    args->solver = (ocp_qp_solver*) c_ptr;
    c_ptr += sizeof(ocp_qp_solver);

    copy_module_pointers_to_args(args->solver, solver);

    assert((size_t)c_ptr % 8 == 0 && "double not 8-byte aligned!");

    args->pcond_args = ocp_qp_partial_condensing_assign_args(dims, c_ptr);
    c_ptr += ocp_qp_partial_condensing_calculate_args_size(dims);

    assert((size_t)c_ptr % 8 == 0 && "double not 8-byte aligned!");

    args->solver_args = args->solver->assign_args(dims, c_ptr);
    c_ptr += args->solver->calculate_args_size(dims);

    assert((char*)raw_memory + ocp_qp_sparse_solver_calculate_args_size(dims, solver) == c_ptr);

    return (void*)args;
}



void ocp_qp_sparse_solver_initialize_default_args(void *args_)
{
    ocp_qp_sparse_solver_args *args = (ocp_qp_sparse_solver_args *)args_;
    ocp_qp_partial_condensing_initialize_default_args(args->pcond_args);
    args->solver->initialize_default_args(args->solver_args);
}



int ocp_qp_sparse_solver_calculate_memory_size(ocp_qp_dims *dims, void *args_)
{
    ocp_qp_sparse_solver_args *args = (ocp_qp_sparse_solver_args *)args_;

    int size = 0;
    size += sizeof(ocp_qp_sparse_solver_memory);

    // set up dimesions of partially condensed qp
    ocp_qp_dims *pcond_dims;
    if (args->pcond_args->N2 < dims->N)
    {
        pcond_dims = args->pcond_args->pcond_dims;
    } else
    {
        pcond_dims = dims;
    }

    if (args->pcond_args->N2 < dims->N) {
        size += ocp_qp_partial_condensing_calculate_memory_size(dims, args->pcond_args);
    }

    size += args->solver->calculate_memory_size(pcond_dims, args->solver_args);

    if (args->pcond_args->N2 < dims->N) {
        size += ocp_qp_in_calculate_size(pcond_dims);
        size += ocp_qp_out_calculate_size(pcond_dims);
    }

    return size;
}



void *ocp_qp_sparse_solver_assign_memory(ocp_qp_dims *dims, void *args_, void *raw_memory)
{
    ocp_qp_sparse_solver_args *args = (ocp_qp_sparse_solver_args *)args_;

    assert(args->pcond_args->N2 > 0 && "N2 must be positive!");
    assert(args->pcond_args->N2 <= dims->N && "N2 cannot be bigger than N!");

    char *c_ptr = (char *)raw_memory;

    // set up dimesions of partially condensed qp
    ocp_qp_dims *pcond_dims;
    if (args->pcond_args->N2 < dims->N)
    {
        pcond_dims = args->pcond_args->pcond_dims;
    } else
    {
        pcond_dims = dims;
    }

    ocp_qp_sparse_solver_memory *mem = (ocp_qp_sparse_solver_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_sparse_solver_memory);

    assert((size_t)c_ptr % 8 == 0 && "double not 8-byte aligned!");

    if (args->pcond_args->N2 < dims->N) {
        mem->pcond_memory = ocp_qp_partial_condensing_assign_memory(dims, args->pcond_args, c_ptr);
        c_ptr += ocp_qp_partial_condensing_calculate_memory_size(dims, args->pcond_args);
    } else
    {
        mem->pcond_memory = NULL;
    }

    assert((size_t)c_ptr % 8 == 0 && "double not 8-byte aligned!");

    mem->solver_memory = args->solver->assign_memory(pcond_dims, args->solver_args, c_ptr);
    c_ptr += args->solver->calculate_memory_size(pcond_dims, args->solver_args);

    if (args->pcond_args->N2 < dims->N) {
        mem->pcond_qp_in = assign_ocp_qp_in(pcond_dims, c_ptr);
        c_ptr += ocp_qp_in_calculate_size(pcond_dims);
    } else
    {
        mem->pcond_qp_in = NULL;
    }

    if (args->pcond_args->N2 < dims->N) {
        mem->pcond_qp_out = assign_ocp_qp_out(pcond_dims, c_ptr);
        c_ptr += ocp_qp_out_calculate_size(pcond_dims);
    } else {
        mem->pcond_qp_out = NULL;
    }

    assert((char *) raw_memory + ocp_qp_sparse_solver_calculate_memory_size(dims, args_) == c_ptr);

    return mem;
}



int ocp_qp_sparse_solver_calculate_workspace_size(ocp_qp_dims *dims, void *args_)
{
    ocp_qp_sparse_solver_args *args = (ocp_qp_sparse_solver_args *)args_;

    int size = sizeof(ocp_qp_sparse_solver_workspace);
    size += ocp_qp_partial_condensing_calculate_workspace_size(dims, args->pcond_args);

    // set up dimesions of partially condensed qp
    ocp_qp_dims *pcond_dims;
    if (args->pcond_args->N2 < dims->N)
    {
        pcond_dims = args->pcond_args->pcond_dims;
    } else
    {
        pcond_dims = dims;
    }

    size += args->solver->calculate_workspace_size(pcond_dims, args->solver_args);

    return size;
}



static void cast_workspace(ocp_qp_dims *dims, ocp_qp_sparse_solver_args *args, ocp_qp_sparse_solver_memory *mem, ocp_qp_sparse_solver_workspace *work)
{
    ocp_qp_dims *pdims = mem->pcond_qp_in->dim;

    char *c_ptr = (char *) work;

    c_ptr += sizeof(ocp_qp_sparse_solver_workspace);

    work->pcond_work = c_ptr;
    c_ptr += ocp_qp_partial_condensing_calculate_workspace_size(dims, args->pcond_args);

    work->solver_work = c_ptr;
    c_ptr += args->solver->calculate_workspace_size(pdims, args->solver_args);
}



int ocp_qp_sparse_solver(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *mem_, void *work_)
{
    ocp_qp_info *info = (ocp_qp_info *)qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer, cond_timer;
    acados_tic(&tot_timer);

    // cast data structures
    ocp_qp_sparse_solver_args *args = (ocp_qp_sparse_solver_args *) args_;
    ocp_qp_sparse_solver_memory *memory = (ocp_qp_sparse_solver_memory *) mem_;
    ocp_qp_sparse_solver_workspace *work = (ocp_qp_sparse_solver_workspace *) work_;

    // cast workspace
    cast_workspace(qp_in->dim, args, memory, work);

    // condense
    acados_tic(&cond_timer);
    if (args->pcond_args->N2 < qp_in->dim->N) {
        ocp_qp_partial_condensing(qp_in, memory->pcond_qp_in, args->pcond_args, memory->pcond_memory, work->pcond_work);
    } else {
        memory->pcond_qp_in = qp_in;
        memory->pcond_qp_out = qp_out;
    }
    info->condensing_time = acados_toc(&cond_timer);

    // solve qp
    int solver_status = args->solver->fun(memory->pcond_qp_in, memory->pcond_qp_out, args->solver_args, memory->solver_memory, work->solver_work);

    // expand
    acados_tic(&cond_timer);
    if (args->pcond_args->N2 < qp_in->dim->N) {
        ocp_qp_partial_expansion(memory->pcond_qp_out, qp_out, args->pcond_args, memory->pcond_memory, work->pcond_work);
    }
    info->condensing_time += acados_toc(&cond_timer);

    info->total_time = acados_toc(&tot_timer);
    info->solve_QP_time = ((ocp_qp_info *)(memory->pcond_qp_out->misc))->solve_QP_time;
    info->interface_time = ((ocp_qp_info *)(memory->pcond_qp_out->misc))->interface_time;
    // return
    return solver_status;
}
