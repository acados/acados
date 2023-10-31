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

// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// acados
#include "acados/sim/sim_common.h"
#include "acados_c/sim_interface.h"
#include "acados/utils/external_function_generic.h"
#include "acados_c/external_function_interface.h"
// example specific
#include "acados_sim_solver_{{ model.name }}.h"
// mex
#include "mex.h"
#include "mex_macros.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    int acados_size, tmp;
    char fun_name[50] = "sim_set";
    char buffer [300]; // for error messages

    /* RHS */

    // C object
    const mxArray *C_sim = prhs[0];
    long long *ptr;
    // solver
    ptr = (long long *) mxGetData( mxGetField( C_sim, 0, "solver" ) );
    sim_solver *solver = (sim_solver *) ptr[0];
    // config
    ptr = (long long *) mxGetData( mxGetField( C_sim, 0, "config" ) );
    sim_config *config = (sim_config *) ptr[0];
    // dims
    ptr = (long long *) mxGetData( mxGetField( C_sim, 0, "dims" ) );
    void *dims = (void *) ptr[0];
    // in
    ptr = (long long *) mxGetData( mxGetField( C_sim, 0, "in" ) );
    sim_in *in = (sim_in *) ptr[0];
    // capsule
    ptr = (long long *) mxGetData( mxGetField( C_sim, 0, "capsule" ) );
    {{ model.name }}_sim_solver_capsule *capsule = ({{ model.name }}_sim_solver_capsule *) ptr[0];

    // field
    char *field = mxArrayToString( prhs[1] );

    // value
    double *value = mxGetPr( prhs[2] );
    int matlab_size = (int) mxGetNumberOfElements( prhs[2] );

    // check dimension, set value
    if (!strcmp(field, "T"))
    {
        {%- if solver_options.integrator_type == "GNSF" %}
            acados_size = 1;
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            sim_in_set(config, dims, in, field, value);
        {% else %}
            MEX_FIELD_NOT_SUPPORTED_FOR_SOLVER(fun_name, field, "irk_gnsf")
        {% endif %}
    }
    else if (!strcmp(field, "x"))
    {
        sim_dims_get(config, dims, "nx", &acados_size);
        MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
        sim_in_set(config, dims, in, field, value);
    }
    else if (!strcmp(field, "u"))
    {
        sim_dims_get(config, dims, "nu", &acados_size);
        MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
        sim_in_set(config, dims, in, field, value);
    }
    else if (!strcmp(field, "p"))
    {
        acados_size = {{ dims.np }};
        MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
        {{ model.name }}_acados_sim_update_params(capsule, value, acados_size);
    }
    else if (!strcmp(field, "xdot"))
    {
        {%- if solver_options.integrator_type == "IRK" or solver_options.integrator_type == "LIFTED_IRK" %}
            sim_dims_get(config, dims, "nx", &acados_size);
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            sim_solver_set(solver, field, value);
        {% else %}
            MEX_FIELD_ONLY_SUPPORTED_FOR_SOLVER(fun_name, field, "irk");
        {% endif %}
    }
    else if (!strcmp(field, "z"))
    {
        {%- if solver_options.integrator_type == "IRK" or solver_options.integrator_type == "LIFTED_IRK" %}
            sim_dims_get(config, dims, "nz", &acados_size);
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            sim_solver_set(solver, field, value);
        {% else %}
            MEX_FIELD_ONLY_SUPPORTED_FOR_SOLVER(fun_name, field, "irk");
        {% endif %}
    }
    else if (!strcmp(field, "phi_guess"))
    {
        {%- if solver_options.integrator_type == "GNSF" %}
            sim_dims_get(config, dims, "nout", &acados_size);
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            sim_solver_set(solver, field, value);
        {% else %}
            MEX_FIELD_ONLY_SUPPORTED_FOR_SOLVER(fun_name, field, "irk_gnsf");
        {% endif %}
    }
    else if (!strcmp(field, "seed_adj"))
    {
        sim_dims_get(config, dims, "nx", &acados_size);
        // TODO(oj): in C, the backward seed is of dimension nx+nu, I think it should only be nx.
        // sim_dims_get(config, dims, "nu", &tmp);
        // acados_size += tmp;
        MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
        sim_in_set(config, dims, in, field, value);
    }
    else
    {
        MEX_FIELD_NOT_SUPPORTED_SUGGEST(fun_name, field, "T, x, u, p, xdot, z, phi_guess, seed_adj");
    }

    return;
}

