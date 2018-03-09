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

#include "acados/sim/sim_common.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/sim/sim_irk_integrator.h"

#include "acados_c/sim_interface.h"

//external
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "acados/utils/mem.h"

sim_solver_config *sim_config_create(sim_solver_plan *plan)
{
    int bytes = sim_solver_config_calculate_size();
    void *ptr = calloc(1, bytes);
    sim_solver_config *solver_config = sim_solver_config_assign(ptr);

    sim_solver_t solver_name = plan->sim_solver;

    // TODO(dimitris): cath error if solver not compiled
    // printf("\n\nSpecified solver interface not compiled with acados!\n\n");
    switch (solver_name)
    {
        case ERK:
            sim_erk_config_initialize_default(solver_config);
            break;
        case LIFTED_IRK:
            sim_lifted_irk_config_initialize_default(solver_config);
            break;
        case IRK:
            sim_irk_config_initialize_default(solver_config);
            break;
    }
    return solver_config;
}



sim_dims *sim_dims_create()
{
    int bytes = sim_dims_calculate_size();

    void *ptr = calloc(1, bytes);

    sim_dims *dims = sim_dims_assign(ptr);

    return dims;
}



sim_in *sim_in_create(sim_solver_config *config, sim_dims *dims)
{
    int bytes = sim_in_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    sim_in *in = sim_in_assign(config, dims, ptr);

    return in;
}



sim_out *sim_out_create(sim_solver_config *config, sim_dims *dims)
{
    int bytes = sim_out_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    sim_out *out = sim_out_assign(config, dims, ptr);

    return out;
}



void *sim_opts_create(sim_solver_config *config, sim_dims *dims)
{
    int bytes = config->opts_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    void *opts = config->opts_assign(config, dims, ptr);

    config->opts_initialize_default(config, dims, opts);

    return opts;
}



int sim_calculate_size(sim_solver_config *config, sim_dims *dims, void *opts_)
{
    int bytes = sizeof(sim_solver);

    bytes += config->memory_calculate_size(config, dims, opts_);
    bytes += config->workspace_calculate_size(config, dims, opts_);

    return bytes;
}



sim_solver *sim_assign(sim_solver_config *config, sim_dims *dims, void *opts_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    sim_solver *solver = (sim_solver *) c_ptr;
    c_ptr += sizeof(sim_solver);

    solver->config = config;
    solver->dims = dims;
    solver->opts = opts_;

    // TODO(dimitris): CHECK ALIGNMENT!

    solver->mem = config->memory_assign(config, dims, opts_, c_ptr);
    c_ptr += config->memory_calculate_size(config, dims, opts_);

    solver->work = (void *) c_ptr;
    c_ptr += config->workspace_calculate_size(config, dims, opts_);

    assert((char*)raw_memory + sim_calculate_size(config, dims, opts_) == c_ptr);

    return solver;
}



sim_solver *sim_create(sim_solver_config *config, sim_dims *dims, void *opts_)
{
    int bytes = sim_calculate_size(config, dims, opts_);

    void *ptr = calloc(1, bytes);

    sim_solver *solver = sim_assign(config, dims, opts_, ptr);

    return solver;
}



int sim_solve(sim_solver *solver, sim_in *qp_in, sim_out *qp_out)
{
    return solver->config->evaluate(solver->config, qp_in, qp_out, solver->opts, solver->mem, solver->work);
}
