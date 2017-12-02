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

// standard
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
// acados
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"


int sim_in_calculate_size(sim_dims *dims)
{
    int size = sizeof(sim_in);

    int nx = dims->nx;
    int nu = dims->nu;

    size += nx * sizeof(double);  // x
    size += nu * sizeof(double);  // u
    size += nx * (nx+nu) * sizeof(double);  // S_forw (max dimension)
    size += (nx + nu) * sizeof(double);  // S_adj

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



sim_in *assign_sim_in(sim_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    sim_in *in = (sim_in *) c_ptr;
    c_ptr += sizeof(sim_in);

    int nx = dims->nx;
    int nu = dims->nu;
    int NF = nx+nu;

    // TODO(dimitris): USE DIMS INSIDE SIM_IN INSTEAD!
    in->nx = nx;
    in->nu = nu;

    align_char_to(8, &c_ptr);

    assign_double(nx, &in->x, &c_ptr);
    assign_double(nu, &in->u, &c_ptr);
    assign_double(nx * NF, &in->S_forw, &c_ptr);
    assign_double(NF, &in->S_adj, &c_ptr);

    assert((char*)raw_memory + sim_in_calculate_size(dims) >= c_ptr);

    return in;
}



// TODO(dimitris): move to create.c
sim_in *create_sim_in(sim_dims *dims)
{
    int bytes = sim_in_calculate_size(dims);

    void *ptr = acados_malloc(bytes, 1);

    sim_in *in = assign_sim_in(dims, ptr);

    return in;
}



int sim_out_calculate_size(sim_dims *dims)
{
    int size = sizeof(sim_out);

    int nx = dims->nx;
    int nu = dims->nu;
    int NF = nx + nu;
    size += sizeof(sim_info);

    size += nx * sizeof(double);  // xn
    size += nx * NF * sizeof(double);  // S_forw
    size += (nx + nu) * sizeof(double);  // S_adj
    size += ((NF + 1) * NF / 2) * sizeof(double);  // S_hess
    size += NF * sizeof(double);  // grad

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



sim_out *assign_sim_out(sim_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    int nx = dims->nx;
    int nu = dims->nu;
    int NF = nx + nu;

    sim_out *out = (sim_out *) c_ptr;
    c_ptr += sizeof(sim_out);

    out->info = (sim_info *)c_ptr;
    c_ptr += sizeof(sim_info);

    align_char_to(8, &c_ptr);

    assign_double(nx, &out->xn, &c_ptr);
    assign_double(nx * NF, &out->S_forw, &c_ptr);
    assign_double(nx + nu, &out->S_adj, &c_ptr);
    assign_double((NF + 1) * NF / 2, &out->S_hess, &c_ptr);
    assign_double(NF, &out->grad, &c_ptr);

    assert((char*)raw_memory + sim_out_calculate_size(dims) >= c_ptr);

    return out;
}



sim_out *create_sim_out(sim_dims *dims)
{
    int bytes = sim_out_calculate_size(dims);

    void *ptr = acados_malloc(bytes, 1);

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
