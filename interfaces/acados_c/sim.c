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

#include "acados_c/sim.h"

//external
#include <stdlib.h>
#include <assert.h>
#include <string.h>
//acados
#include <acados/utils/mem.h>
//acados_c
#include "acados_c/sim/sim_erk_integrator.h"
#include "acados_c/sim/sim_lifted_irk_integrator.h"



void sim_copy_dims(sim_dims *dest, sim_dims *src)
{
    dest->num_stages = src->num_stages;

    dest->nx = src->nx;
    
    dest->nu = src->nu;
}



sim_dims *create_sim_dims()
{
    int bytes = sim_dims_calculate_size();

    void *ptr = malloc(bytes);

    sim_dims *dims = assign_sim_dims(ptr);

    return dims;
}



sim_in *create_sim_in(sim_dims *dims)
{
    int bytes = sim_in_calculate_size(dims);

    void *ptr = acados_malloc(bytes, 1);

    sim_in *in = assign_sim_in(dims, ptr);

    return in;
}



sim_out *create_sim_out(sim_dims *dims)
{
    int bytes = sim_out_calculate_size(dims);

    void *ptr = malloc(bytes);

    sim_out *out = assign_sim_out(dims, ptr);

    return out;
}



int sim_calculate_args_size(sim_solver_plan *plan, sim_dims *dims)
{
    sim_solver_fcn_ptrs fcn_ptrs;

    set_sim_solver_fcn_ptrs(plan, &fcn_ptrs);

    int size = fcn_ptrs.calculate_args_size(dims);

    return size;
}



void *sim_assign_args(sim_solver_plan *plan, sim_dims *dims, void *raw_memory)
{
    sim_solver_fcn_ptrs fcn_ptrs;

    set_sim_solver_fcn_ptrs(plan, &fcn_ptrs);

    void *args = fcn_ptrs.assign_args(dims, raw_memory);

    fcn_ptrs.initialize_default_args(dims, args);

    return args;
}



void *sim_create_args(sim_solver_plan *plan, sim_dims *dims)
{
    int bytes = sim_calculate_args_size(plan, dims);

    void *ptr = malloc(bytes);

    void *args = sim_assign_args(plan, dims, ptr);

    return args;
}



void *sim_copy_args(sim_solver_plan *plan, sim_dims *dims, void *raw_memory, void *source)
{
    sim_solver_fcn_ptrs fcn_ptrs;

    set_sim_solver_fcn_ptrs(plan, &fcn_ptrs);

    void *args = fcn_ptrs.copy_args(dims, raw_memory, source);

    return args;
}



int sim_calculate_size(sim_solver_plan *plan, sim_dims *dims, void *args_)
{
    sim_solver_fcn_ptrs fcn_ptrs;

    set_sim_solver_fcn_ptrs(plan, &fcn_ptrs);

    int bytes = 0;

    bytes += sizeof(sim_solver);

    bytes += sizeof(sim_solver_fcn_ptrs);

    bytes += sim_dims_calculate_size();

    bytes += fcn_ptrs.calculate_args_size(dims);

    bytes += fcn_ptrs.calculate_memory_size(dims, args_);

    bytes += fcn_ptrs.calculate_workspace_size(dims, args_);

    return bytes;
}



sim_solver *sim_assign(sim_solver_plan *plan, sim_dims *dims, void *args_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    sim_rk_opts *args = (sim_rk_opts *)args_;

    sim_solver *solver = (sim_solver *) c_ptr;
    c_ptr += sizeof(sim_solver);

    solver->fcn_ptrs = (sim_solver_fcn_ptrs *) c_ptr;
    c_ptr += sizeof(sim_solver_fcn_ptrs);
    set_sim_solver_fcn_ptrs(plan, solver->fcn_ptrs);

    solver->dims = assign_sim_dims(c_ptr);
    c_ptr += sim_dims_calculate_size();
    sim_copy_dims(solver->dims, dims);

    solver->args = solver->fcn_ptrs->copy_args(dims, c_ptr, args_);
    c_ptr += solver->fcn_ptrs->calculate_args_size(dims);

    solver->mem = solver->fcn_ptrs->assign_memory(dims, args_, c_ptr);
    c_ptr += solver->fcn_ptrs->calculate_memory_size(dims, args_);

    solver-> work = (void *) c_ptr;
    c_ptr += solver->fcn_ptrs->calculate_workspace_size(dims, args_);

    assert((char*)raw_memory + sim_calculate_size(plan, dims, args_) == c_ptr);

    return solver;
}



sim_solver *sim_create(sim_solver_plan *plan, sim_dims *dims, void *args_)
{
    int bytes = sim_calculate_size(plan, dims, args_);

    void *ptr = malloc(bytes);

    sim_solver *solver = sim_assign(plan, dims, args_, ptr);

    return solver;
}



int sim_solve(sim_solver *solver, sim_in *qp_in, sim_out *qp_out)
{
    return solver->fcn_ptrs->fun(qp_in, qp_out, solver->args, solver->mem, solver->work);
}



int set_sim_solver_fcn_ptrs(sim_solver_plan *plan, sim_solver_fcn_ptrs *fcn_ptrs)
{
    int return_value = ACADOS_SUCCESS;
    sim_solver_t solver_name = plan->sim_solver;

    switch (solver_name)
    {
        case ERK:
            fcn_ptrs->fun = &sim_erk;
            fcn_ptrs->calculate_args_size = &sim_erk_opts_calculate_size;
            fcn_ptrs->assign_args = &sim_erk_assign_opts;
            fcn_ptrs->initialize_default_args = &sim_erk_initialize_default_args;
            fcn_ptrs->calculate_memory_size = &sim_erk_calculate_memory_size;
            fcn_ptrs->assign_memory = &sim_erk_assign_memory;
            fcn_ptrs->calculate_workspace_size = &sim_erk_calculate_workspace_size;
            break;
        case LIFTED_IRK:
            fcn_ptrs->fun = &sim_lifted_irk;
            fcn_ptrs->calculate_args_size = &sim_lifted_irk_opts_calculate_size;
            fcn_ptrs->assign_args = &sim_lifted_irk_assign_opts;
            fcn_ptrs->initialize_default_args = &sim_lifted_irk_initialize_default_args;
            fcn_ptrs->calculate_memory_size = &sim_lifted_irk_calculate_memory_size;
            fcn_ptrs->assign_memory = &sim_lifted_irk_assign_memory;
            fcn_ptrs->calculate_workspace_size = &sim_lifted_irk_calculate_workspace_size;
            break;
        default:
            return_value = ACADOS_FAILURE;
    }

    return return_value;
}