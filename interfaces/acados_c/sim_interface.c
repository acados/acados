/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
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

sim_config *sim_config_create(sim_solver_plan_t plan)
{
    /* calculate_size & malloc & assign */

    acados_size_t bytes = sim_config_calculate_size();
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
        case INVALID_SIM_SOLVER:
            printf("\nerror: sim_config_create: forgot to initialize plan->sim_solver\n");
            exit(1);
            break;
        default:
            printf("\nerror: sim_config_create: unsupported plan->sim_solver\n");
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
//    acados_size_t bytes = config->dims_calculate_size(config_);
    acados_size_t bytes = config->dims_calculate_size();

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


void sim_dims_get_from_attr(sim_config *config, void *dims, const char *field, int *dims_out)
{
    // vectors first
    dims_out[1] = 0;

    if (!strcmp(field, "x") || !strcmp(field, "xdot"))
    {
        sim_dims_get(config, dims, "nx", &dims_out[0]);
    }
    else if (!strcmp(field, "z"))
    {
        sim_dims_get(config, dims, "nz", &dims_out[0]);
    }
    else if (!strcmp(field, "u"))
    {
        sim_dims_get(config, dims, "nu", &dims_out[0]);
    }
    else if (!strcmp(field, "T") || !strcmp(field, "t0"))
    {
        dims_out[0] = 1;
    }
    else if (!strcmp(field, "S_adj"))
    {
        int tmp;
        sim_dims_get(config, dims, "nu", &tmp);
        sim_dims_get(config, dims, "nx", &dims_out[0]);
        dims_out[0] += tmp;
    }
    else if (!strcmp(field, "S_algebraic"))
    {
        int tmp;
        sim_dims_get(config, dims, "nz", &dims_out[0]);
        sim_dims_get(config, dims, "nx", &dims_out[1]);
        sim_dims_get(config, dims, "nu", &tmp);
        dims_out[1] += tmp;
    }
    else if (!strcmp(field, "S_forw"))
    {
        sim_dims_get(config, dims, "nx", &dims_out[0]);
        sim_dims_get(config, dims, "nu", &dims_out[1]);
        dims_out[1] += dims_out[0];
    }
    else if (!strcmp(field, "S_hess"))
    {
        sim_dims_get(config, dims, "nx", &dims_out[0]);
        sim_dims_get(config, dims, "nu", &dims_out[1]);
        dims_out[1] += dims_out[0];
        dims_out[0] = dims_out[1];
    }
    else if (!strcmp(field, "Sx"))
    {
        sim_dims_get(config, dims, "nx", &dims_out[0]);
        sim_dims_get(config, dims, "nx", &dims_out[1]);
    }
    else if (!strcmp(field, "Su"))
    {
        sim_dims_get(config, dims, "nx", &dims_out[0]);
        sim_dims_get(config, dims, "nu", &dims_out[1]);
    }
    else if (!strcmp(field, "seed_adj"))
    {
        sim_dims_get(config, dims, "nx", &dims_out[0]);
    }
    else
    {
        printf("\nerror: sim_dims_get_from_attr: field %s not available\n", field);
        exit(1);
    }

    return;
}

/************************************************
* in
************************************************/

sim_in *sim_in_create(sim_config *config, void *dims)
{
    acados_size_t bytes = sim_in_calculate_size(config, dims);

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
    acados_size_t bytes = sim_out_calculate_size(config, dims);

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
    acados_size_t bytes = config->opts_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    void *opts = config->opts_assign(config, dims, ptr);

    config->opts_initialize_default(config, dims, opts);

    return opts;
}



void sim_opts_destroy(void *opts)
{
    free(opts);
}


void sim_opts_set(sim_config *config, void *opts, const char *field,
                           void *value)
{
    config->opts_set(config, opts, field, value);
}

void sim_opts_get(sim_config *config, void *opts, const char *field, void *value)
{
    sim_opts *opts_ = (sim_opts *) opts;
    sim_opts_get_(config, opts_, field, value);
}

/************************************************
* solver
************************************************/

acados_size_t sim_calculate_size(sim_config *config, void *dims, void *opts_)
{
    acados_size_t bytes = sizeof(sim_solver);

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
    acados_size_t bytes = sim_calculate_size(config, dims, opts_);

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


int sim_solver_set(sim_solver *solver, const char *field, void *value)
{
    return solver->config->memory_set(solver->config, solver->dims, solver->mem, field, value);
}
