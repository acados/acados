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
#include "acados_c/sim_interface.h"
// example specific
#include "acados_sim_solver_{{ name }}.h"

// mex
#include "mex.h"
#include "mex_macros.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    // sizeof(long long) == sizeof(void *) = 64 !!!
    long long *l_ptr;
    char fun_name[50] = "sim_create";
    int status = 0;

    // create sim solver
    {{ name }}_sim_solver_capsule *acados_sim_capsule = {{ name }}_acados_sim_solver_create_capsule();
    status = {{ name }}_acados_sim_create(acados_sim_capsule);
    if (status)
    {
        mexPrintf("{{ name }}_acados_create() returned status %d.\n", status);
    }
    mexPrintf("{{ name }}_acados_create() -> success!\n");

    /* RHS */
    // no input params

    /* LHS */
    #define FIELDS_SIM 8

    // field names of output struct
    char *fieldnames[FIELDS_SIM];
    fieldnames[0] = (char*)mxMalloc(50);
    fieldnames[1] = (char*)mxMalloc(50);
    fieldnames[2] = (char*)mxMalloc(50);
    fieldnames[3] = (char*)mxMalloc(50);
    fieldnames[4] = (char*)mxMalloc(50);
    fieldnames[5] = (char*)mxMalloc(50);
    fieldnames[6] = (char*)mxMalloc(50);
    fieldnames[7] = (char*)mxMalloc(50);

    memcpy(fieldnames[0],"config",sizeof("config"));
    memcpy(fieldnames[1],"dims",sizeof("dims"));
    memcpy(fieldnames[2],"opts",sizeof("opts"));
    memcpy(fieldnames[3],"in",sizeof("in"));
    memcpy(fieldnames[4],"out",sizeof("out"));
    memcpy(fieldnames[5],"solver",sizeof("solver"));
    memcpy(fieldnames[6],"capsule",sizeof("capsule"));
    memcpy(fieldnames[7],"mem",sizeof("mem"));

    // create output struct
    plhs[0] = mxCreateStructMatrix(1, 1, FIELDS_SIM, (const char **) fieldnames);

    mxFree( fieldnames[0] );
    mxFree( fieldnames[1] );
    mxFree( fieldnames[2] );
    mxFree( fieldnames[3] );
    mxFree( fieldnames[4] );
    mxFree( fieldnames[5] );
    mxFree( fieldnames[6] );
    mxFree( fieldnames[7] );


    /* populate output struct */
    // config
    mxArray *config_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(config_mat);
    sim_config * config = {{ name }}_acados_get_sim_config(acados_sim_capsule);
    l_ptr[0] = (long long) config;
    mxSetField(plhs[0], 0, "config", config_mat);

    // dims
    mxArray *dims_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(dims_mat);
    void * dims = {{ name }}_acados_get_sim_dims(acados_sim_capsule);
    l_ptr[0] = (long long) dims;
    mxSetField(plhs[0], 0, "dims", dims_mat);

    // opts
    mxArray *opts_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(opts_mat);
    sim_opts * opts = {{ name }}_acados_get_sim_opts(acados_sim_capsule);
    l_ptr[0] = (long long) opts;
    mxSetField(plhs[0], 0, "opts", opts_mat);

    // in
    mxArray *in_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(in_mat);
    sim_in * in = {{ name }}_acados_get_sim_in(acados_sim_capsule);
    l_ptr[0] = (long long) in;
    mxSetField(plhs[0], 0, "in", in_mat);

    // out
    mxArray *out_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(out_mat);
    sim_out * out = {{ name }}_acados_get_sim_out(acados_sim_capsule);
    l_ptr[0] = (long long) out;
    mxSetField(plhs[0], 0, "out", out_mat);

    // solver
    mxArray *solver_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(solver_mat);
    sim_solver * solver = {{ name }}_acados_get_sim_solver(acados_sim_capsule);
    l_ptr[0] = (long long) solver;
    mxSetField(plhs[0], 0, "solver", solver_mat);

    // capsule
    mxArray *capsule_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(capsule_mat);
    l_ptr[0] = (long long) acados_sim_capsule;
    mxSetField(plhs[0], 0, "capsule", capsule_mat);

    // mem
    mxArray *mem_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(mem_mat);
    void * mem = {{ name }}_acados_get_sim_mem(acados_sim_capsule);
    l_ptr[0] = (long long) mem;
    mxSetField(plhs[0], 0, "mem", mem_mat);

    return;

}
