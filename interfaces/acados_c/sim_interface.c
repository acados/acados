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
#include "acados/sim/sim_gnsf.h"
#include "acados/sim/sim_irk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"

#include "acados_c/sim_interface.h"

// external
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "acados/utils/mem.h"



/************************************************
* config
************************************************/

sim_config *sim_config_create(sim_solver_plan plan)
{
    /* calculate_size & malloc & assign */

    int bytes = sim_config_calculate_size();
    void *ptr = calloc(1, bytes);
    sim_config *solver_config = sim_config_assign(ptr);

    /* initialize config according plan */

    sim_solver_t solver_name = plan.sim_solver;

    switch (solver_name)
    {
        case ERK:
            sim_erk_config_initialize_default(solver_config);
            break;
        case IRK:
            sim_irk_config_initialize_default(solver_config);
            break;
        case GNSF:
            sim_gnsf_config_initialize_default(solver_config);
            break;
        case LIFTED_IRK:
            sim_lifted_irk_config_initialize_default(solver_config);
            break;
        default:
            printf("\n\nSpecified integrator not available in acados C interface!\n\n");
            exit(1);
    }
    return solver_config;
}



void sim_config_destroy(void *config)
{
    free(config);
}



/************************************************
* dims
************************************************/

void *sim_dims_create(void *config_)
{
    sim_config *config = (sim_config *) config_;
//    int bytes = config->dims_calculate_size(config_);
    int bytes = config->dims_calculate_size();

    void *ptr = calloc(1, bytes);

    void *dims = config->dims_assign(config_, ptr);

    return dims;
}



void sim_dims_destroy(void *dims)
{
    free(dims);
}



void sim_dims_set(sim_config *config, void *dims, const char *field, const int* value)
{
    config->dims_set(config, dims, field, value);
}



void sim_dims_get(sim_config *config, void *dims, const char *field, int* value)
{
    config->dims_get(config, dims, field, value);
}



/************************************************
* in
************************************************/

sim_in *sim_in_create(sim_config *config, void *dims)
{
    int bytes = sim_in_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    sim_in *in = sim_in_assign(config, dims, ptr);

    return in;
}



void sim_in_destroy(void *in)
{
    free(in);
}



int sim_in_set(void *config_, void *dims_, sim_in *in, const char *field, void *value)
{
    return sim_in_set_(config_, dims_, in, field, value);
}



/************************************************
* out
************************************************/


sim_out *sim_out_create(sim_config *config, void *dims)
{
    int bytes = sim_out_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    sim_out *out = sim_out_assign(config, dims, ptr);

    return out;
}



void sim_out_destroy(void *out)
{
    free(out);
}



int sim_out_get(void *config, void *dims, sim_out *out, const char *field, void *value)
{
    return sim_out_get_(config, dims, out, field, value);
}


/************************************************
* options
************************************************/

void *sim_opts_create(sim_config *config, void *dims)
{
    int bytes = config->opts_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    void *opts = config->opts_assign(config, dims, ptr);

    config->opts_initialize_default(config, dims, opts);

    return opts;
}



void sim_opts_destroy(void *opts)
{
    free(opts);
}


int sim_opts_set(sim_config *config, void *opts, const char *field,
                           void *value)
{
    return config->opts_set(config, opts, field, value);
}


/************************************************
* solver
************************************************/

int sim_calculate_size(sim_config *config, void *dims, void *opts_)
{
    int bytes = sizeof(sim_solver);

    bytes += config->memory_calculate_size(config, dims, opts_);
    bytes += config->workspace_calculate_size(config, dims, opts_);

    return bytes;
}



sim_solver *sim_assign(sim_config *config, void *dims, void *opts_, void *raw_memory)
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

    assert((char *) raw_memory + sim_calculate_size(config, dims, opts_) == c_ptr);

    return solver;
}



sim_solver *sim_solver_create(sim_config *config, void *dims, void *opts_)
{
    // update Butcher tableau (needed if the user changed ns)
    config->opts_update(config, dims, opts_);
    int bytes = sim_calculate_size(config, dims, opts_);

    void *ptr = calloc(1, bytes);

    sim_solver *solver = sim_assign(config, dims, opts_, ptr);

    return solver;
}



void sim_solver_destroy(void *solver)
{
    free(solver);
}



int sim_solve(sim_solver *solver, sim_in *in, sim_out *out)
{
    int status = solver->config->evaluate(solver->config, in, out, solver->opts, solver->mem,
                                    solver->work);
    return status;
}

int sim_precompute(sim_solver *solver, sim_in *in, sim_out *out)
{
    return solver->config->precompute(solver->config, in, out, solver->opts, solver->mem,
                                    solver->work);
}
