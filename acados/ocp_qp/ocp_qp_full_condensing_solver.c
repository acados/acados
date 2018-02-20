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
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"



int ocp_qp_full_condensing_solver_opts_calculate_size(void *config_, ocp_qp_dims *dims)
{
	ocp_qp_xcond_solver_config *config = config_;
	dense_qp_solver_config *qp_solver = config->qp_solver;

    int size = 0;
    size += sizeof(ocp_qp_full_condensing_solver_args);

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    size += ocp_qp_full_condensing_calculate_args_size(dims);
    size += qp_solver->opts_calculate_size(config_, &ddims);

    return size;
}



void *ocp_qp_full_condensing_solver_opts_assign(void *config_, ocp_qp_dims *dims, void *raw_memory)
{
	ocp_qp_xcond_solver_config *config = config_;
	dense_qp_solver_config *qp_solver = config->qp_solver;

    char *c_ptr = (char *) raw_memory;

    ocp_qp_full_condensing_solver_args *args = (ocp_qp_full_condensing_solver_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_full_condensing_solver_args);

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    args->cond_args = ocp_qp_full_condensing_assign_args(dims, c_ptr);
    c_ptr += ocp_qp_full_condensing_calculate_args_size(dims);

	align_char_to(8, &c_ptr);

    args->solver_args = qp_solver->opts_assign(qp_solver, &ddims, c_ptr);
    c_ptr += qp_solver->opts_calculate_size(qp_solver, &ddims);

    assert((char*)raw_memory + ocp_qp_full_condensing_solver_opts_calculate_size(config_, dims) == c_ptr);

    return (void*)args;
}



void ocp_qp_full_condensing_solver_opts_initialize_default(void *config_, void *args_)
{
	ocp_qp_xcond_solver_config *config = config_;
	dense_qp_solver_config *qp_solver = config->qp_solver;

	// full cond solver
    ocp_qp_full_condensing_solver_args *args = (ocp_qp_full_condensing_solver_args *)args_;
	// full condensing
    ocp_qp_full_condensing_initialize_default_args(args->cond_args);
	// qp solver
    qp_solver->opts_initialize_default(qp_solver, args->solver_args);
}



int ocp_qp_full_condensing_solver_memory_calculate_size(void *config_, ocp_qp_dims *dims, void *args_)
{
	ocp_qp_xcond_solver_config *config = config_;
	dense_qp_solver_config *qp_solver = config->qp_solver;

    ocp_qp_full_condensing_solver_args *args = (ocp_qp_full_condensing_solver_args *)args_;

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    int size = 0;
    size += sizeof(ocp_qp_full_condensing_solver_memory);

    size += ocp_qp_full_condensing_calculate_memory_size(dims, args->cond_args);
    size += qp_solver->memory_calculate_size(qp_solver, &ddims, args->solver_args);

    size += dense_qp_in_calculate_size(&ddims);
    size += dense_qp_out_calculate_size(&ddims);

    return size;
}



void *ocp_qp_full_condensing_solver_memory_assign(void *config_, ocp_qp_dims *dims, void *args_, void *raw_memory)
{
	ocp_qp_xcond_solver_config *config = config_;
	dense_qp_solver_config *qp_solver = config->qp_solver;

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

    mem->solver_memory = qp_solver->memory_assign(qp_solver, &ddims, args->solver_args, c_ptr);
    c_ptr += qp_solver->memory_calculate_size(qp_solver, &ddims, args->solver_args);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    mem->qpd_in = assign_dense_qp_in(&ddims, c_ptr);
    c_ptr += dense_qp_in_calculate_size(&ddims);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    mem->qpd_out = assign_dense_qp_out(&ddims, c_ptr);
    c_ptr += dense_qp_out_calculate_size(&ddims);

    assert((char *) raw_memory + ocp_qp_full_condensing_solver_memory_calculate_size(config_, dims, args_) == c_ptr);

    return mem;
}



int ocp_qp_full_condensing_solver_workspace_calculate_size(void *config_, ocp_qp_dims *dims, void *args_)
{
	ocp_qp_xcond_solver_config *config = config_;
	dense_qp_solver_config *qp_solver = config->qp_solver;

    ocp_qp_full_condensing_solver_args *args = (ocp_qp_full_condensing_solver_args *)args_;

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    int size = sizeof(ocp_qp_full_condensing_solver_workspace);
    size += ocp_qp_full_condensing_calculate_workspace_size(dims, args->cond_args);
    size += qp_solver->workspace_calculate_size(qp_solver, &ddims, args->solver_args);

    return size;
}



static void cast_workspace(void *config_, ocp_qp_dims *dims, ocp_qp_full_condensing_solver_args *args, ocp_qp_full_condensing_solver_memory *mem, ocp_qp_full_condensing_solver_workspace *work)
{
	ocp_qp_xcond_solver_config *config = config_;
	dense_qp_solver_config *qp_solver = config->qp_solver;

    dense_qp_dims *ddims = mem->qpd_in->dim;

    char *c_ptr = (char *) work;

    c_ptr += sizeof(ocp_qp_full_condensing_solver_workspace);

    work->cond_work = c_ptr;
    c_ptr += ocp_qp_full_condensing_calculate_workspace_size(dims, args->cond_args);

    work->solver_workspace = c_ptr;
    c_ptr += qp_solver->workspace_calculate_size(qp_solver, ddims, args->solver_args);
}



int ocp_qp_full_condensing_solver(void *config_, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *mem_, void *work_)
{
	ocp_qp_xcond_solver_config *config = config_;
	dense_qp_solver_config *qp_solver = config->qp_solver;

    ocp_qp_info *info = (ocp_qp_info *)qp_out->misc;
    acados_timer tot_timer, cond_timer;
    acados_tic(&tot_timer);

    // cast data structures
    ocp_qp_full_condensing_solver_args *args = (ocp_qp_full_condensing_solver_args *) args_;
    ocp_qp_full_condensing_solver_memory *memory = (ocp_qp_full_condensing_solver_memory *) mem_;
    ocp_qp_full_condensing_solver_workspace *work = (ocp_qp_full_condensing_solver_workspace *) work_;

    // cast workspace
    cast_workspace(config_, qp_in->dim, args, memory, work);

    // condense
    acados_tic(&cond_timer);
    ocp_qp_full_condensing(qp_in, memory->qpd_in, args->cond_args, memory->cond_memory, work->cond_work);
    info->condensing_time = acados_toc(&cond_timer);

    // solve qp
    int solver_status = qp_solver->evaluate(qp_solver, memory->qpd_in, memory->qpd_out, args->solver_args, memory->solver_memory, work->solver_workspace);

    // expand
    acados_tic(&cond_timer);
    ocp_qp_full_expansion(memory->qpd_out, qp_out, args->cond_args, memory->cond_memory, work->cond_work);
    info->condensing_time += acados_toc(&cond_timer);

    info->total_time = acados_toc(&tot_timer);
    info->solve_QP_time = ((dense_qp_info *)(memory->qpd_out->misc))->solve_QP_time;
    info->interface_time = ((dense_qp_info *)(memory->qpd_out->misc))->interface_time;
    info->num_iter = ((dense_qp_info *)(memory->qpd_out->misc))->num_iter;

    return solver_status;
}



void ocp_qp_full_condensing_solver_config_initialize_default(void *config_)
{

	ocp_qp_xcond_solver_config *config = config_;

	config->opts_calculate_size = &ocp_qp_full_condensing_solver_opts_calculate_size;
	config->opts_assign = &ocp_qp_full_condensing_solver_opts_assign;
	config->opts_initialize_default = &ocp_qp_full_condensing_solver_opts_initialize_default;
	config->memory_calculate_size = &ocp_qp_full_condensing_solver_memory_calculate_size;
	config->memory_assign = &ocp_qp_full_condensing_solver_memory_assign;
	config->workspace_calculate_size = &ocp_qp_full_condensing_solver_workspace_calculate_size;
	config->evaluate = &ocp_qp_full_condensing_solver;

	return;

}
