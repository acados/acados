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

// mex
#include "mex.h"
#include "mex_macros.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    int acados_size, tmp;
    long long *ptr;
    char fun_name[50] = "sim_set";
    char buffer [300]; // for error messages

    /* RHS */

    // model_struct
    char *ext_fun_type;
    const mxArray *matlab_model = prhs[0];
    if (mxGetField( matlab_model, 0, "ext_fun_type" )!=NULL)
        ext_fun_type = mxArrayToString( mxGetField( matlab_model, 0, "ext_fun_type" ) );

    const mxArray *C_sim = prhs[2];
    const mxArray *C_ext_fun_pointers = prhs[3];

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
    // method
    int method = mxGetScalar( mxGetField( C_sim, 0, "method" ) );

    // field
    char *field = mxArrayToString( prhs[4] );
    double *value = mxGetPr( prhs[5] );

    int matlab_size = (int) mxGetNumberOfElements( prhs[5] );

    // check dimension, set value
    if (!strcmp(field, "T"))
    {
        if (method != GNSF)
        {
            acados_size = 1;
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            sim_in_set(config, dims, in, field, value);
        }
        else
        {
            MEX_FIELD_NOT_SUPPORTED_FOR_SOLVER(fun_name, field, "irk_gnsf")
        }
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
        // loop over number of external functions;
        int struct_size = mxGetNumberOfFields( C_ext_fun_pointers );
        for (int ii=0; ii<struct_size; ii++)
        {
//            printf("\n%s\n", mxGetFieldNameByNumber( C_ext_fun_pointers, ii) );
            ptr = (long long *) mxGetData( mxGetFieldByNumber( C_ext_fun_pointers, 0, ii ) );
            // external function param casadi
            if (!strcmp(ext_fun_type, "casadi"))
            {
                if (ptr[0]!=0)
                {
                    external_function_param_casadi *ext_fun_ptr = (external_function_param_casadi *) ptr[0];
                    acados_size = ext_fun_ptr->np;
                    MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);

                    ext_fun_ptr->set_param(ext_fun_ptr, value);
                }
            }
            // external function param generic
            else if (!strcmp(ext_fun_type, "generic"))
            {
                if (ptr[0]!=0)
                {
                    external_function_param_generic *ext_fun_ptr = (external_function_param_generic *) ptr[0];
                    acados_size = ext_fun_ptr->np;
                    MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);

                    ext_fun_ptr->set_param(ext_fun_ptr, value);
                }
            }
            else
            {
                MEX_FIELD_VALUE_NOT_SUPPORTED_SUGGEST(fun_name, "ext_fun_type", ext_fun_type, "casadi, generic");
            }
        }
    }
    else if (!strcmp(field, "xdot"))
    {
        if (method == IRK || method == LIFTED_IRK)
        {
            sim_dims_get(config, dims, "nx", &acados_size);
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            sim_solver_set(solver, field, value);
        }
        else
        {
            MEX_FIELD_ONLY_SUPPORTED_FOR_SOLVER(fun_name, field, "irk");
        }
    }
    else if (!strcmp(field, "z"))
    {
        if (method == IRK || method == LIFTED_IRK)
        {
            sim_dims_get(config, dims, "nz", &acados_size);
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            sim_solver_set(solver, field, value);
        }
        else
        {
            MEX_FIELD_ONLY_SUPPORTED_FOR_SOLVER(fun_name, field, "irk");
        }
    }
    else if (!strcmp(field, "phi_guess"))
    {
        if (method == GNSF)
        {
            sim_dims_get(config, dims, "nout", &acados_size);
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            sim_solver_set(solver, field, value);
        }
        else
        {
            MEX_FIELD_ONLY_SUPPORTED_FOR_SOLVER(fun_name, field, "irk_gnsf");
        }
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

