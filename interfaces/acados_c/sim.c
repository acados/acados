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



int set_sim_solver_fun_ptrs(sim_solver_t sim_solver_name, sim_solver *sim_solver)
{
    int return_value = ACADOS_SUCCESS;

    switch (sim_solver_name)
    {
        case ERK:
            sim_solver->fun = &sim_erk;
            sim_solver->calculate_args_size = &sim_erk_opts_calculate_size;
            sim_solver->assign_args = &sim_erk_assign_opts;
            sim_solver->initialize_default_args = &sim_erk_initialize_default_args;
            sim_solver->calculate_memory_size = &sim_erk_calculate_memory_size;
            sim_solver->assign_memory = &sim_erk_assign_memory;
            sim_solver->calculate_workspace_size = &sim_erk_calculate_workspace_size;
            break;
        case LIFTED_IRK:
            sim_solver->fun = &sim_lifted_irk;
            sim_solver->calculate_args_size = &sim_lifted_irk_opts_calculate_size;
            sim_solver->assign_args = &sim_lifted_irk_assign_opts;
            sim_solver->initialize_default_args = &sim_lifted_irk_initialize_default_args;
            sim_solver->calculate_memory_size = &sim_lifted_irk_calculate_memory_size;
            sim_solver->assign_memory = &sim_lifted_irk_assign_memory;
            sim_solver->calculate_workspace_size = &sim_lifted_irk_calculate_workspace_size;
            break;
        default:
            return_value = ACADOS_FAILURE;
    }
    return return_value;
}