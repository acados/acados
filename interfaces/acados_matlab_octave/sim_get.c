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
// mex
#include "mex.h"
#include "mex_macros.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    long long *ptr;
    char fun_name[20] = "sim_get";
    char buffer [200]; // for error messages

    /* RHS */

    // C_sim

    // config
    ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "config" ) );
    sim_config *config = (sim_config *) ptr[0];
    // dims
    ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "dims" ) );
    void *dims = (void *) ptr[0];
    // out
    ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "out" ) );
    sim_out *out = (sim_out *) ptr[0];

    // field
    char *field = mxArrayToString( prhs[1] );
//    mexPrintf("\nin sim_get: field"\n%s\n", field);

    int nx; sim_dims_get(config, dims, "nx", &nx);
    int nu; sim_dims_get(config, dims, "nu", &nu);
    int nz; sim_dims_get(config, dims, "nz", &nz);


    if (!strcmp(field, "xn") || !strcmp(field, "x"))
    {
        plhs[0] = mxCreateNumericMatrix(nx, 1, mxDOUBLE_CLASS, mxREAL);
        double *xn = mxGetPr( plhs[0] );
        sim_out_get(config, dims, out, "xn", xn);
    }
    else if (!strcmp(field, "zn") || !strcmp(field, "z"))
    {
        plhs[0] = mxCreateNumericMatrix(nz, 1, mxDOUBLE_CLASS, mxREAL);
        double *zn = mxGetPr( plhs[0] );
        sim_out_get(config, dims, out, "zn", zn);
    }
    else if (!strcmp(field, "S_forw"))
    {
        plhs[0] = mxCreateNumericMatrix(nx, nu+nx, mxDOUBLE_CLASS, mxREAL);
        double *S_forw = mxGetPr( plhs[0] );
        sim_out_get(config, dims, out, "S_forw", S_forw);
    }
    else if (!strcmp(field, "Sx"))
    {
        plhs[0] = mxCreateNumericMatrix(nx, nx, mxDOUBLE_CLASS, mxREAL);
        double *Sx = mxGetPr( plhs[0] );
        sim_out_get(config, dims, out, "Sx", Sx);
    }
    else if (!strcmp(field, "Su"))
    {
        plhs[0] = mxCreateNumericMatrix(nx, nu, mxDOUBLE_CLASS, mxREAL);
        double *Su = mxGetPr( plhs[0] );
        sim_out_get(config, dims, out, "Su", Su);
    }
    else if (!strcmp(field, "S_adj"))
    {
        plhs[0] = mxCreateNumericMatrix(nx+nu, 1, mxDOUBLE_CLASS, mxREAL);
        double *S_adj = mxGetPr( plhs[0] );
        sim_out_get(config, dims, out, "S_adj", S_adj);
    }
    else if (!strcmp(field, "S_hess"))
    {
        plhs[0] = mxCreateNumericMatrix(nx+nu, nx+nu, mxDOUBLE_CLASS, mxREAL);
        double *S_hess = mxGetPr( plhs[0] );
        sim_out_get(config, dims, out, "S_hess", S_hess);
    }
    else if (!strcmp(field, "S_algebraic"))
    {
        plhs[0] = mxCreateNumericMatrix(nz, nx+nu, mxDOUBLE_CLASS, mxREAL);
        double *S_algebraic = mxGetPr( plhs[0] );
        sim_out_get(config, dims, out, "S_algebraic", S_algebraic);
    }
    else if (!strcmp(field, "time_tot"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *time_tot = mxGetPr( plhs[0] );
        sim_out_get(config, dims, out, "time_tot", time_tot);
    }
    else if (!strcmp(field, "time_la"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *time_la = mxGetPr( plhs[0] );
        sim_out_get(config, dims, out, "time_la", time_la);
    }
    else if (!strcmp(field, "time_ad"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *time_ad = mxGetPr( plhs[0] );
        sim_out_get(config, dims, out, "time_ad", time_ad);
    }
    else
    {
        MEX_FIELD_NOT_SUPPORTED_SUGGEST(fun_name, field,
             "xn, x, zn, z, S_forw, Sx, Su, S_hess, S_algebraic, time_tot, time_la, time_ad");
    }

    return;

}


