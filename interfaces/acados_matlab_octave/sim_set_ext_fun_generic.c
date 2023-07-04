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



// external functions for the model
int FUN_NAME(void **, void **, void *);



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // sizeof(long long) == sizeof(void *) = 64 !!!
    long long *ptr;

    /* RHS */
    const mxArray *C_sim = prhs[0];
    const mxArray *matlab_model = prhs[2];

    // config
    ptr = (long long *) mxGetData( mxGetField( C_sim, 0, "config" ) );
    sim_config *config = (sim_config *) ptr[0];
    // dims
    ptr = (long long *) mxGetData( mxGetField( C_sim, 0, "dims" ) );
    void *dims = (void *) ptr[0];
    // in
    ptr = (long long *) mxGetData( mxGetField( C_sim, 0, "in" ) );
    sim_in *in = (sim_in *) ptr[0];

    // model
    int np = 0;
    if (mxGetField( matlab_model, 0, "dim_np" )!=NULL)
    {
        np = mxGetScalar( mxGetField( matlab_model, 0, "dim_np" ) );
    }

    /* LHS */
    /* copy existing fields */
    plhs[0] = mxDuplicateArray(prhs[1]);

    /* populate new fields */
    external_function_param_generic *ext_fun_ptr;

    ext_fun_ptr = (external_function_param_generic *) malloc(1*sizeof(external_function_param_generic));
    external_function_param_generic_set_fun(ext_fun_ptr, &FUN_NAME);
    external_function_param_generic_create(ext_fun_ptr, np);
    sim_in_set(config, dims, in, STR(SET_FIELD), ext_fun_ptr);
    // populate output struct
    ptr = mxGetData(mxGetField(plhs[0], 0, STR(MEX_FIELD)));
    ptr[0] = (long long) ext_fun_ptr;

    return;

}



