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



int sim_calculate_args_size(sim_solver_fcn_ptrs *fcn_ptrs, sim_dims *dims)
{
    return fcn_ptrs->calculate_args_size(dims, fcn_ptrs->submodules);
}



void *sim_assign_args(sim_solver_fcn_ptrs *fcn_ptrs, sim_dims *dims, void *raw_memory)
{
    void *args = fcn_ptrs->assign_args(dims, &fcn_ptrs->submodules, raw_memory);

    fcn_ptrs->initialize_default_args(dims, args);

    return args;
}



void *sim_create_args(sim_solver_fcn_ptrs *fcn_ptrs, sim_dims *dims)
{
    int bytes = sim_calculate_args_size(fcn_ptrs, dims);

    void *ptr = malloc(bytes);

    void *args = sim_assign_args(fcn_ptrs, dims, ptr);

    return args;
}



void *sim_copy_args(sim_solver_fcn_ptrs *fcn_ptrs, sim_dims *dims, void *raw_memory, void *source)
{
    return fcn_ptrs->copy_args(dims, raw_memory, source);
}



int sim_calculate_size(sim_solver_fcn_ptrs *fcn_ptrs, sim_dims *dims, void *args_)
{
    int bytes = 0;

    bytes += sizeof(sim_solver);

    bytes += sizeof(sim_solver_fcn_ptrs);

    bytes += sim_dims_calculate_size();

    bytes += sim_calculate_args_size(fcn_ptrs, dims);

    bytes += fcn_ptrs->calculate_memory_size(dims, args_);

    bytes += fcn_ptrs->calculate_workspace_size(dims, args_);

    return bytes;
}



sim_solver *sim_assign(sim_solver_fcn_ptrs *fcn_ptrs, sim_dims *dims, void *args_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    sim_solver *solver = (sim_solver *) c_ptr;
    c_ptr += sizeof(sim_solver);

    solver->fcn_ptrs = (sim_solver_fcn_ptrs *) c_ptr;
    c_ptr += sizeof(sim_solver_fcn_ptrs);

    solver->dims = assign_sim_dims(c_ptr);
    c_ptr += sim_dims_calculate_size();

    solver->args = sim_copy_args(fcn_ptrs, dims, c_ptr, args_);
    c_ptr += sim_calculate_args_size(fcn_ptrs, dims);

    solver->mem = fcn_ptrs->assign_memory(dims, args_, c_ptr);
    c_ptr += fcn_ptrs->calculate_memory_size(dims, args_);

    solver->work = (void *) c_ptr;
    c_ptr += fcn_ptrs->calculate_workspace_size(dims, args_);

    assert((char*)raw_memory + sim_calculate_size(fcn_ptrs, dims, args_) == c_ptr);

    *solver->fcn_ptrs = *fcn_ptrs;
    solver->fcn_ptrs->submodules = NULL;

    sim_copy_dims(solver->dims, dims);

    return solver;
}



sim_solver *sim_create(sim_solver_fcn_ptrs *fcn_ptrs, sim_dims *dims, void *args_)
{
    int bytes = sim_calculate_size(fcn_ptrs, dims, args_);

    void *ptr = malloc(bytes);

    sim_solver *solver = sim_assign(fcn_ptrs, dims, args_, ptr);

    return solver;
}



int sim_solve(sim_solver *solver, sim_in *qp_in, sim_out *qp_out)
{
    return solver->fcn_ptrs->fun(qp_in, qp_out, solver->args, solver->mem, solver->work);
}



int sim_calculate_submodules_size(sim_solver_config *config, sim_dims *dims)
{
    sim_solver_t solver_name = config->sim_solver;

    int size;

    switch (solver_name)
    {
        case ERK:
            size = sim_erk_integrator_calculate_submodules_size(config, dims);
            break;
        case LIFTED_IRK:
            size = sim_lifted_irk_integrator_calculate_submodules_size(config, dims);
            break;
        default:
            size = 0;
    }

    return size;
}



void *sim_assign_submodules(sim_solver_config *config, sim_dims *dims, void *raw_memory)
{
    sim_solver_t solver_name = config->sim_solver;

    void *submodules;

    switch (solver_name)
    {
        case ERK:
            submodules = sim_erk_integrator_assign_submodules(config, dims, raw_memory);
            break;
        case LIFTED_IRK:
            submodules = sim_lifted_irk_integrator_assign_submodules(config, dims, raw_memory);
            break;
        default:
            submodules = NULL;
    }

    return submodules;
}



int calculate_sim_solver_fcn_ptrs_size(sim_solver_config *config, sim_dims *dims)
{
    int size = sizeof(sim_solver_fcn_ptrs);

    size += sim_calculate_submodules_size(config, dims);

    return size;
}



void *assign_sim_solver_fcn_ptrs(sim_solver_config *config, sim_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *)raw_memory;

    sim_solver_fcn_ptrs *fcn_ptrs = (sim_solver_fcn_ptrs *)c_ptr;
    c_ptr += sizeof(sim_solver_fcn_ptrs);

    set_sim_solver_fcn_ptrs(config, fcn_ptrs);

    fcn_ptrs->submodules = sim_assign_submodules(config, dims, c_ptr);
    c_ptr += sim_calculate_submodules_size(config, dims);

    assert((char*)raw_memory + calculate_sim_solver_fcn_ptrs_size(config, dims) == c_ptr);

    return (void *)fcn_ptrs;
}



void *create_sim_solver_fcn_ptrs(sim_solver_config *config, sim_dims *dims)
{
    int bytes = calculate_sim_solver_fcn_ptrs_size(config, dims);

    void *ptr = malloc(bytes);

    sim_solver_fcn_ptrs *fcn_ptrs = assign_sim_solver_fcn_ptrs(config, dims, ptr);

    return fcn_ptrs;
}



int set_sim_solver_fcn_ptrs(sim_solver_config *config, sim_solver_fcn_ptrs *fcn_ptrs)
{
    int return_value = ACADOS_SUCCESS;
    sim_solver_t solver_name = config->sim_solver;

    switch (solver_name)
    {
        case ERK:
            fcn_ptrs->fun = &sim_erk_integrator;
            fcn_ptrs->calculate_args_size = &sim_erk_integrator_calculate_args_size;
            fcn_ptrs->assign_args = &sim_erk_integrator_assign_args;
            fcn_ptrs->copy_args = &sim_erk_integrator_copy_args;
            fcn_ptrs->initialize_default_args = &sim_erk_integrator_initialize_default_args;
            fcn_ptrs->calculate_memory_size = &sim_erk_integrator_calculate_memory_size;
            fcn_ptrs->assign_memory = &sim_erk_integrator_assign_memory;
            fcn_ptrs->calculate_workspace_size = &sim_erk_integrator_calculate_workspace_size;
            break;
        case LIFTED_IRK:
            fcn_ptrs->fun = &sim_lifted_irk_integrator;
            fcn_ptrs->calculate_args_size = &sim_lifted_irk_integrator_calculate_args_size;
            fcn_ptrs->assign_args = &sim_lifted_irk_integrator_assign_args;
            fcn_ptrs->copy_args = &sim_lifted_irk_integrator_copy_args;
            fcn_ptrs->initialize_default_args = &sim_lifted_irk_integrator_initialize_default_args;
            fcn_ptrs->calculate_memory_size = &sim_lifted_irk_integrator_calculate_memory_size;
            fcn_ptrs->assign_memory = &sim_lifted_irk_integrator_assign_memory;
            fcn_ptrs->calculate_workspace_size = &sim_lifted_irk_integrator_calculate_workspace_size;
            break;
        default:
            return_value = ACADOS_FAILURE;
    }

    return return_value;
}
