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


int ocp_qp_sparse_solver_calculate_args_size(ocp_qp_dims *dims, void *submodules_)
{
    ocp_qp_sparse_solver_submodules *submodules = (ocp_qp_sparse_solver_submodules *) submodules_;

    int size = 0;
    size += sizeof(ocp_qp_sparse_solver_args);
    size += sizeof(ocp_qp_solver_fcn_ptrs);

    size += ocp_qp_partial_condensing_calculate_args_size(dims, NULL);
    size += submodules->solver->calculate_args_size(dims, submodules->solver->submodules);

    return size;
}



void *ocp_qp_sparse_solver_assign_args(ocp_qp_dims *dims, void **submodules_, void *raw_memory)
{
    ocp_qp_sparse_solver_submodules *submodules = (ocp_qp_sparse_solver_submodules *) *submodules_;

    char *c_ptr = (char *) raw_memory;

    ocp_qp_sparse_solver_args *args = (ocp_qp_sparse_solver_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_sparse_solver_args);

    args->submodules.solver = (ocp_qp_solver_fcn_ptrs*) c_ptr;
    c_ptr += sizeof(ocp_qp_solver_fcn_ptrs);

    assert((size_t)c_ptr % 8 == 0 && "double not 8-byte aligned!");

    void *pcond_submodules = NULL;
    args->pcond_args = ocp_qp_partial_condensing_assign_args(dims, &pcond_submodules, c_ptr);
    c_ptr += ocp_qp_partial_condensing_calculate_args_size(dims, NULL);

    assert((size_t)c_ptr % 8 == 0 && "double not 8-byte aligned!");

    // QUESTION(nielsvd): why dims and not pcond_dims?
    void *solver_submodules = submodules->solver->submodules;
    args->solver_args = submodules->solver->assign_args(dims, &solver_submodules, c_ptr);
    c_ptr += submodules->solver->calculate_args_size(dims, submodules->solver->submodules);

    assert((char*)raw_memory + ocp_qp_sparse_solver_calculate_args_size(dims, *submodules_) == c_ptr);

    // Update submodules' fcn_ptrs
    *(args->submodules.solver) = *(submodules->solver);
    // Update submodules' submodules pointer
    args->submodules.solver->submodules = solver_submodules;

    // Update submodules pointer
    *submodules_ = (void *)&(args->submodules);

    return (void*)args;
}



void *ocp_qp_sparse_solver_copy_args(ocp_qp_dims *dims, void *raw_memory, void *source_)
{
    ocp_qp_sparse_solver_args *source = (ocp_qp_sparse_solver_args *) source_;
    ocp_qp_sparse_solver_args *dest;

    ocp_qp_sparse_solver_submodules *submodules = &source->submodules;

    dest = ocp_qp_sparse_solver_assign_args(dims, (void **)&submodules, raw_memory);

    ocp_qp_partial_condensing_copy_args(dims, dest->pcond_args, source->pcond_args);

    source->submodules.solver->copy_args(dims, dest->solver_args, source->solver_args);

    return (void *) dest;
}



void ocp_qp_sparse_solver_initialize_default_args(void *args_)
{
    ocp_qp_sparse_solver_args *args = (ocp_qp_sparse_solver_args *)args_;
    ocp_qp_partial_condensing_initialize_default_args(args->pcond_args);
    args->submodules.solver->initialize_default_args(args->solver_args);
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
    }
	else
    {
        pcond_dims = dims;
    }

    if (args->pcond_args->N2 < dims->N) {
        size += ocp_qp_partial_condensing_calculate_memory_size(dims, args->pcond_args);
    }

    size += args->submodules.solver->calculate_memory_size(pcond_dims, args->solver_args);

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
    }
	else
    {
        pcond_dims = dims;
    }

    ocp_qp_sparse_solver_memory *mem = (ocp_qp_sparse_solver_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_sparse_solver_memory);

    assert((size_t)c_ptr % 8 == 0 && "double not 8-byte aligned!");

    if (args->pcond_args->N2 < dims->N) {
        mem->pcond_memory = ocp_qp_partial_condensing_assign_memory(dims, args->pcond_args, c_ptr);
        c_ptr += ocp_qp_partial_condensing_calculate_memory_size(dims, args->pcond_args);
    }
	else
    {
        mem->pcond_memory = NULL;
    }

    assert((size_t)c_ptr % 8 == 0 && "double not 8-byte aligned!");

    mem->solver_memory = args->submodules.solver->assign_memory(pcond_dims, args->solver_args, c_ptr);
    c_ptr += args->submodules.solver->calculate_memory_size(pcond_dims, args->solver_args);

    if (args->pcond_args->N2 < dims->N) {
        mem->pcond_qp_in = assign_ocp_qp_in(pcond_dims, c_ptr);
        c_ptr += ocp_qp_in_calculate_size(pcond_dims);
    }
	else
    {
        mem->pcond_qp_in = NULL;
    }

    if (args->pcond_args->N2 < dims->N) {
        mem->pcond_qp_out = assign_ocp_qp_out(pcond_dims, c_ptr);
        c_ptr += ocp_qp_out_calculate_size(pcond_dims);
    }
	else
	{
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
    }
	else
    {
        pcond_dims = dims;
    }

    size += args->submodules.solver->calculate_workspace_size(pcond_dims, args->solver_args);

    return size;
}



static void cast_workspace(ocp_qp_dims *dims, ocp_qp_sparse_solver_args *args, ocp_qp_sparse_solver_memory *mem, ocp_qp_sparse_solver_workspace *work)
{
    // set up dimesions of partially condensed qp
    ocp_qp_dims *pcond_dims;
    if (args->pcond_args->N2 < dims->N)
    {
        pcond_dims = args->pcond_args->pcond_dims;
    }
	else
    {
        pcond_dims = dims;
    }

    char *c_ptr = (char *) work;

    c_ptr += sizeof(ocp_qp_sparse_solver_workspace);

    work->pcond_work = c_ptr;
    c_ptr += ocp_qp_partial_condensing_calculate_workspace_size(dims, args->pcond_args);

    work->solver_work = c_ptr;
    c_ptr += args->submodules.solver->calculate_workspace_size(pcond_dims, args->solver_args);
}



int ocp_qp_sparse_solver(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *mem_, void *work_)
{
    ocp_qp_info *info = (ocp_qp_info *)qp_out->misc;
    acados_timer tot_timer, cond_timer;
    acados_tic(&tot_timer);

    // cast data structures
    ocp_qp_sparse_solver_args *args = (ocp_qp_sparse_solver_args *) args_;
    ocp_qp_sparse_solver_memory *memory = (ocp_qp_sparse_solver_memory *) mem_;
    ocp_qp_sparse_solver_workspace *work = (ocp_qp_sparse_solver_workspace *) work_;

    // cast workspace
    cast_workspace(qp_in->dim, args, memory, work);

    if (args->pcond_args->N2 < qp_in->dim->N) {  // condensing
        acados_tic(&cond_timer);
        ocp_qp_partial_condensing(qp_in, memory->pcond_qp_in, args->pcond_args, memory->pcond_memory, work->pcond_work);
        info->condensing_time = acados_toc(&cond_timer);
    } else {
        memory->pcond_qp_in = qp_in;
        memory->pcond_qp_out = qp_out;
        info->condensing_time = 0;
    }

    // solve qp
    int solver_status = args->submodules.solver->fun(memory->pcond_qp_in, memory->pcond_qp_out, args->solver_args, memory->solver_memory, work->solver_work);

    // expand
    if (args->pcond_args->N2 < qp_in->dim->N) {
        acados_tic(&cond_timer);
        ocp_qp_partial_expansion(memory->pcond_qp_out, qp_out, args->pcond_args, memory->pcond_memory, work->pcond_work);
        info->condensing_time += acados_toc(&cond_timer);
    }

    info->total_time = acados_toc(&tot_timer);
    info->solve_QP_time = ((ocp_qp_info *)(memory->pcond_qp_out->misc))->solve_QP_time;
    info->interface_time = ((ocp_qp_info *)(memory->pcond_qp_out->misc))->interface_time;
    info->num_iter = ((ocp_qp_info *)(memory->pcond_qp_out->misc))->num_iter;

    return solver_status;
}
