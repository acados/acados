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
#include "acados_c/sim/sim_irk_integrator.h"
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

    sim_dims *dims = sim_dims_assign(ptr);

    return dims;
}



sim_in *create_sim_in(sim_dims *dims, sim_solver_config *config)
{
    int bytes = sim_in_calculate_size(dims, config);

    void *ptr = acados_malloc(bytes, 1);

    sim_in *in = sim_in_assign(dims, ptr, config);

    return in;
}



sim_out *create_sim_out(sim_dims *dims)
{
    int bytes = sim_out_calculate_size(dims);

    void *ptr = malloc(bytes);

    sim_out *out = sim_out_assign(dims, ptr);

    return out;
}



int sim_calculate_args_size(sim_solver_plan *plan, sim_dims *dims)
{
    sim_solver_config fcn_ptrs = {
        /*.fun = */ NULL,
        /*.calculate_args_size = */ NULL,
        /*.assign_args = */ NULL,
        /*.initialize_default_args = */ NULL,
        /*.calculate_memory_size = */ NULL,
        /*.assign_memory = */ NULL,
        /*.calculate_workspace_size = */ NULL};

    set_sim_solver_fcn_ptrs(plan, &fcn_ptrs);

    int size = fcn_ptrs.opts_calculate_size(dims);

    return size;
}



void *sim_assign_args(sim_solver_plan *plan, sim_dims *dims, void *raw_memory)
{
    sim_solver_config fcn_ptrs = {
        /*.fun = */ NULL,
        /*.calculate_args_size = */ NULL,
        /*.assign_args = */ NULL,
        /*.initialize_default_args = */ NULL,
        /*.calculate_memory_size = */ NULL,
        /*.assign_memory = */ NULL,
        /*.calculate_workspace_size = */ NULL};

    set_sim_solver_fcn_ptrs(plan, &fcn_ptrs);

    void *args = fcn_ptrs.opts_assign(dims, raw_memory);

    fcn_ptrs.opts_initialize_default(dims, args);

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
    sim_solver_t solver_name = plan->sim_solver;

    void *args = NULL;

    switch (solver_name)
    {
        case ERK:
            args = sim_erk_copy_opts(dims, raw_memory, source);
            break;
        case IRK:
            // TODO(dimitris): NOT IMPLEMENTED YET
            exit(1);
            // args = sim_irk_copy_opts(dims, raw_memory, source);
            break;
        case LIFTED_IRK:
            args = sim_lifted_irk_copy_opts(dims, raw_memory, source);
            break;
//		TODO
//        case IRK:
//            args = sim_lifted_irk_copy_opts(dims, raw_memory, source);
    }

    return args;
}



int sim_calculate_size(sim_solver_plan *plan, sim_dims *dims, void *args_)
{
    sim_solver_config fcn_ptrs = {
        /*.fun = */ NULL,
        /*.calculate_args_size = */ NULL,
        /*.assign_args = */ NULL,
        /*.initialize_default_args = */ NULL,
        /*.calculate_memory_size = */ NULL,
        /*.assign_memory = */ NULL,
        /*.calculate_workspace_size = */ NULL};

    set_sim_solver_fcn_ptrs(plan, &fcn_ptrs);

    int bytes = 0;

    bytes += sizeof(sim_solver);

    bytes += sizeof(sim_solver_config);

    bytes += sim_dims_calculate_size();

    bytes += fcn_ptrs.opts_calculate_size(dims);

    bytes += fcn_ptrs.memory_calculate_size(dims, args_);

    bytes += fcn_ptrs.workspace_calculate_size(dims, args_);

    return bytes;
}



sim_solver *sim_assign(sim_solver_plan *plan, sim_dims *dims, void *args_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    sim_solver *solver = (sim_solver *) c_ptr;
    c_ptr += sizeof(sim_solver);

    solver->fcn_ptrs = (sim_solver_config *) c_ptr;
    c_ptr += sizeof(sim_solver_config);
    set_sim_solver_fcn_ptrs(plan, solver->fcn_ptrs);

    solver->dims = sim_dims_assign(c_ptr);
    c_ptr += sim_dims_calculate_size();
    sim_copy_dims(solver->dims, dims);

    solver->args = sim_copy_args(plan, dims, c_ptr, args_);
    c_ptr += solver->fcn_ptrs->opts_calculate_size(dims);

    solver->mem = solver->fcn_ptrs->memory_assign(dims, args_, c_ptr);
    c_ptr += solver->fcn_ptrs->memory_calculate_size(dims, args_);

    solver-> work = (void *) c_ptr;
    c_ptr += solver->fcn_ptrs->workspace_calculate_size(dims, args_);

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



int set_sim_solver_fcn_ptrs(sim_solver_plan *plan, sim_solver_config *fcn_ptrs)
{
    int return_value = ACADOS_SUCCESS;
    sim_solver_t solver_name = plan->sim_solver;

    switch (solver_name)
    {
        case ERK:
			sim_erk_config_initialize_default(fcn_ptrs);
            break;
        case LIFTED_IRK:
			sim_lifted_irk_config_initialize_default(fcn_ptrs);
            break;
        case IRK:
			sim_irk_config_initialize_default(fcn_ptrs);
            break;
        default:
            return_value = ACADOS_FAILURE;
    }

    return return_value;
}
