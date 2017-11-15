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
#include "acados/ocp_qp/ocp_qp_partial_condensing_solver.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/dense_qp/dense_qp_qpoases.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"


int ocp_qp_partial_condensing_solver_calculate_args_size(ocp_qp_dims *dims, ocp_qp_solver *solver) {

    int size = 0;
    size += sizeof(ocp_qp_partial_condensing_solver_args);

    size += ocp_qp_partial_condensing_calculate_args_size(dims);
    size += solver->calculate_args_size(dims);

    return size;
}



void *ocp_qp_partial_condensing_solver_assign_args(ocp_qp_dims *dims, ocp_qp_solver *solver, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_partial_condensing_solver_args *args = (ocp_qp_partial_condensing_solver_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_partial_condensing_solver_args);

    args->pcond_args = ocp_qp_partial_condensing_assign_args(dims, c_ptr);
    c_ptr += ocp_qp_partial_condensing_calculate_args_size(dims);

    args->solver = solver;

    args->solver_args = solver->assign_args(dims, c_ptr);
    c_ptr += solver->calculate_args_size(dims);

    assert((char*)raw_memory + ocp_qp_partial_condensing_solver_calculate_args_size(dims, solver) >= c_ptr);

    return (void*)args;
}



void ocp_qp_partial_condensing_solver_initialize_default_args(void *args_)
{
    ocp_qp_partial_condensing_solver_args *args = (ocp_qp_partial_condensing_solver_args *)args_;
    ocp_qp_partial_condensing_initialize_default_args(args->pcond_args);
    args->solver->initialize_default_args(args->solver_args);
}



int ocp_qp_partial_condensing_solver_calculate_memory_size(ocp_qp_dims *dims, void *args_)
{
    ocp_qp_partial_condensing_solver_args *args = (ocp_qp_partial_condensing_solver_args *)args_;

    int size = 0;
    size += sizeof(ocp_qp_partial_condensing_memory);

    size += ocp_qp_partial_condensing_calculate_memory_size(dims, args->pcond_args);
    size += args->solver->calculate_memory_size(args->pcond_args->pcond_dims, args->solver_args);

    size += ocp_qp_in_calculate_size(args->pcond_args->pcond_dims);
    size += ocp_qp_out_calculate_size(args->pcond_args->pcond_dims);

    make_int_multiple_of(8, &size);
    size += 4 * 8;

    return size;
}



void *ocp_qp_partial_condensing_solver_assign_memory(ocp_qp_dims *dims, void *args_, void *raw_memory)
{
    ocp_qp_partial_condensing_solver_args *args = (ocp_qp_partial_condensing_solver_args *)args_;

    char *c_ptr = (char *)raw_memory;

    ocp_qp_partial_condensing_solver_memory *mem = (ocp_qp_partial_condensing_solver_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_partial_condensing_solver_memory);

    align_char_to(8, &c_ptr);
    mem->pcond_memory = ocp_qp_partial_condensing_assign_memory(dims, args->pcond_args, c_ptr);
    c_ptr += ocp_qp_partial_condensing_calculate_memory_size(dims, args->pcond_args);

    align_char_to(8, &c_ptr);
    mem->solver_memory = args->solver->assign_memory(args->pcond_args->pcond_dims, args->solver_args, c_ptr);
    c_ptr += args->solver->calculate_memory_size(args->pcond_args->pcond_dims, args->solver_args);

    align_char_to(8, &c_ptr);
    mem->pcond_qp_in = assign_ocp_qp_in(args->pcond_args->pcond_dims, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(args->pcond_args->pcond_dims);

    align_char_to(8, &c_ptr);
    mem->pcond_qp_out = assign_ocp_qp_out(args->pcond_args->pcond_dims, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(args->pcond_args->pcond_dims);

    assert((char *) raw_memory + ocp_qp_partial_condensing_solver_calculate_memory_size(dims, args_) >= c_ptr);

    return mem;
}



int ocp_qp_partial_condensing_solver(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *mem_) {

    ocp_qp_partial_condensing_solver_args *args = (ocp_qp_partial_condensing_solver_args *) args_;
    ocp_qp_partial_condensing_solver_memory *memory = (ocp_qp_partial_condensing_solver_memory *) mem_;

    // condense
    ocp_qp_partial_condensing(qp_in, memory->pcond_qp_in, args->pcond_args, memory->pcond_memory);

    // solve qp
    int solver_status = args->solver->fun(memory->pcond_qp_in, memory->pcond_qp_out, args->solver_args, memory->solver_memory);

    // expand
    ocp_qp_partial_expansion(memory->pcond_qp_out, qp_out, args->pcond_args, memory->pcond_memory);

    // return
    return solver_status;
}
