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
//acados
#include <acados/sim/sim_erk_integrator.h>
#include <acados/sim/sim_lifted_irk_integrator.h>
#include <acados/utils/mem.h>



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
    return 0;
}



void *sim_assign_args(sim_solver_plan *plan, sim_dims *dims, void *raw_memory)
{
    return NULL;
}



void *sim_create_args(sim_solver_plan *plan, sim_dims *dims)
{
    return NULL;
}



int sim_calculate_size(sim_solver_plan *plan, sim_dims *dims, void *args_)
{
    return 0;
}



sim_solver *sim_assign(sim_solver_plan *plan, sim_dims *dims, void *args_, void *raw_memory)
{
    return NULL;
}



sim_solver *sim_create(sim_solver_plan *plan, sim_dims *dims, void *args_)
{
    return NULL;
}



int sim_solve(sim_solver *solver, sim_in *qp_in, sim_out *qp_out)
{
    return 0;
}



int set_sim_solver_fun_ptrs(sim_solver_plan *plan, sim_solver_fcn_ptrs *fcn_ptrs)
{
    int return_value = ACADOS_SUCCESS;
    sim_solver_t sim_solver_name = plan->sim_solver;

    switch (sim_solver_name)
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