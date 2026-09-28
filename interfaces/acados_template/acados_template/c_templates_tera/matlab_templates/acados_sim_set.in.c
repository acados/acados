/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
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
#include "acados_sim_solver_{{ name }}.h"
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
    {{ name }}_sim_solver_capsule *capsule = ({{ name }}_sim_solver_capsule *) ptr[0];

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
        {{ name }}_acados_sim_update_params(capsule, value, acados_size);
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

