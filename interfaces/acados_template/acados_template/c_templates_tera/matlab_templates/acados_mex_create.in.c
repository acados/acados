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


// standard
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// acados
#include "acados/utils/print.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_solver_{{ model.name }}.h"

// mex
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    long long *l_ptr;
    int status = 0;

    // create solver
    {{ model.name }}_solver_capsule *acados_ocp_capsule = {{ model.name }}_acados_create_capsule();

    status = {{ model.name }}_acados_create(acados_ocp_capsule);

    if (status)
    {
        mexPrintf("{{ model.name }}_acados_create() returned status %d.\n", status);
    }
    mexPrintf("{{ model.name }}_acados_create() -> success!\n");

    // get pointers to nlp solver related objects
    ocp_nlp_plan_t *nlp_plan = {{ model.name }}_acados_get_nlp_plan(acados_ocp_capsule);
    ocp_nlp_config *nlp_config = {{ model.name }}_acados_get_nlp_config(acados_ocp_capsule);
    ocp_nlp_dims *nlp_dims = {{ model.name }}_acados_get_nlp_dims(acados_ocp_capsule);
    ocp_nlp_in *nlp_in = {{ model.name }}_acados_get_nlp_in(acados_ocp_capsule);
    ocp_nlp_out *nlp_out = {{ model.name }}_acados_get_nlp_out(acados_ocp_capsule);
    ocp_nlp_solver *nlp_solver = {{ model.name }}_acados_get_nlp_solver(acados_ocp_capsule);
    void *nlp_opts = {{ model.name }}_acados_get_nlp_opts(acados_ocp_capsule);

    // mexPrintf("acados: got pointer to objectes!\n");

    // field names of output struct
    #define FIELDS_OCP 9
    #define FIELDS_EXT_FUN 25
    #define MAX_FIELDS 25
    char *fieldnames[MAX_FIELDS];

    for (int i = 0; i < MAX_FIELDS; i++)
    {
        fieldnames[i] = (char*) mxMalloc(50);
    }

    memcpy(fieldnames[0],"config",sizeof("config"));
    memcpy(fieldnames[1],"dims",sizeof("dims"));
    memcpy(fieldnames[2],"opts",sizeof("opts"));
    memcpy(fieldnames[3],"in",sizeof("in"));
    memcpy(fieldnames[4],"out",sizeof("out"));
    memcpy(fieldnames[5],"solver",sizeof("solver"));
    memcpy(fieldnames[6],"sens_out",sizeof("sens_out"));
    memcpy(fieldnames[7],"plan",sizeof("plan"));
    memcpy(fieldnames[8],"capsule",sizeof("capsule"));

    // create output struct - C_ocp
    plhs[0] = mxCreateStructMatrix(1, 1, 9, (const char **) fieldnames);

    // MEX: config, dims, opts, in, out, solver, sens_out, plan
    // plan
    mxArray *plan_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(plan_mat);
    l_ptr[0] = (long long) nlp_plan;
    mxSetField(plhs[0], 0, "plan", plan_mat);

    // config
    mxArray *config_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(config_mat);
    l_ptr[0] = (long long) nlp_config;
    mxSetField(plhs[0], 0, "config", config_mat);

    // dims
    mxArray *dims_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(dims_mat);
    l_ptr[0] = (long long) nlp_dims;
    mxSetField(plhs[0], 0, "dims", dims_mat);

    // opts
    mxArray *opts_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(opts_mat);
    l_ptr[0] = (long long) nlp_opts;
    mxSetField(plhs[0], 0, "opts", opts_mat);

    // in
    mxArray *in_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(in_mat);
    l_ptr[0] = (long long) nlp_in;
    mxSetField(plhs[0], 0, "in", in_mat);

    // out
    mxArray *out_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(out_mat);
    l_ptr[0] = (long long) nlp_out;
    mxSetField(plhs[0], 0, "out", out_mat);

    // solver
    mxArray *solver_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(solver_mat);
    l_ptr[0] = (long long) nlp_solver;
    mxSetField(plhs[0], 0, "solver", solver_mat);

    // TODO: sens_out not actually implemented in templates..
    // sens_out
    mxArray *sens_out_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(sens_out_mat);
    l_ptr[0] = (long long) 1;
    mxSetField(plhs[0], 0, "sens_out", sens_out_mat);

    // capsule
    mxArray *capsule_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(capsule_mat);
    l_ptr[0] = (long long) acados_ocp_capsule;
    mxSetField(plhs[0], 0, "capsule", capsule_mat);

    return;
}
