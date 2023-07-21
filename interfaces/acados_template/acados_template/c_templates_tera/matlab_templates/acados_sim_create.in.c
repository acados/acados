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
#include "acados_c/sim_interface.h"
// example specific
#include "acados_sim_solver_{{ model.name }}.h"

// mex
#include "mex.h"
#include "mex_macros.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    // sizeof(long long) == sizeof(void *) = 64 !!!
    int nx, nu, ii;

    long long *l_ptr;
    char *c_ptr;
    char fun_name[50] = "sim_create";
    char buffer [300]; // for error messages
    int status = 0;

    // create sim solver
    {{ model.name }}_sim_solver_capsule *acados_sim_capsule = {{ model.name }}_acados_sim_solver_create_capsule();
    status = {{ model.name }}_acados_sim_create(acados_sim_capsule);
    if (status)
    {
        mexPrintf("{{ model.name }}_acados_create() returned status %d.\n", status);
    }
    mexPrintf("{{ model.name }}_acados_create() -> success!\n");

    /* RHS */
    // no input params

    /* LHS */
    #define FIELDS_SIM 7

    // field names of output struct
    char *fieldnames[FIELDS_SIM];
    fieldnames[0] = (char*)mxMalloc(50);
    fieldnames[1] = (char*)mxMalloc(50);
    fieldnames[2] = (char*)mxMalloc(50);
    fieldnames[3] = (char*)mxMalloc(50);
    fieldnames[4] = (char*)mxMalloc(50);
    fieldnames[5] = (char*)mxMalloc(50);
    fieldnames[6] = (char*)mxMalloc(50);

    memcpy(fieldnames[0],"config",sizeof("config"));
    memcpy(fieldnames[1],"dims",sizeof("dims"));
    memcpy(fieldnames[2],"opts",sizeof("opts"));
    memcpy(fieldnames[3],"in",sizeof("in"));
    memcpy(fieldnames[4],"out",sizeof("out"));
    memcpy(fieldnames[5],"solver",sizeof("solver"));
    memcpy(fieldnames[6],"capsule",sizeof("capsule"));

    // create output struct
    plhs[0] = mxCreateStructMatrix(1, 1, FIELDS_SIM, (const char **) fieldnames);

    mxFree( fieldnames[0] );
    mxFree( fieldnames[1] );
    mxFree( fieldnames[2] );
    mxFree( fieldnames[3] );
    mxFree( fieldnames[4] );
    mxFree( fieldnames[5] );
    mxFree( fieldnames[6] );


    /* populate output struct */
    // config
    mxArray *config_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(config_mat);
    sim_config * config = {{ model.name }}_acados_get_sim_config(acados_sim_capsule);
    l_ptr[0] = (long long) config;
    mxSetField(plhs[0], 0, "config", config_mat);

    // dims
    mxArray *dims_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(dims_mat);
    void * dims = {{ model.name }}_acados_get_sim_dims(acados_sim_capsule);
    l_ptr[0] = (long long) dims;
    mxSetField(plhs[0], 0, "dims", dims_mat);

    // opts
    mxArray *opts_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(opts_mat);
    sim_opts * opts = {{ model.name }}_acados_get_sim_opts(acados_sim_capsule);
    l_ptr[0] = (long long) opts;
    mxSetField(plhs[0], 0, "opts", opts_mat);

    // in
    mxArray *in_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(in_mat);
    sim_in * in = {{ model.name }}_acados_get_sim_in(acados_sim_capsule);
    l_ptr[0] = (long long) in;
    mxSetField(plhs[0], 0, "in", in_mat);

    // out
    mxArray *out_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(out_mat);
    sim_out * out = {{ model.name }}_acados_get_sim_out(acados_sim_capsule);
    l_ptr[0] = (long long) out;
    mxSetField(plhs[0], 0, "out", out_mat);

    // solver
    mxArray *solver_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(solver_mat);
    sim_solver * solver = {{ model.name }}_acados_get_sim_solver(acados_sim_capsule);
    l_ptr[0] = (long long) solver;
    mxSetField(plhs[0], 0, "solver", solver_mat);

    // capsule
    mxArray *capsule_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(capsule_mat);
    l_ptr[0] = (long long) acados_sim_capsule;
    mxSetField(plhs[0], 0, "capsule", capsule_mat);

    return;

}
