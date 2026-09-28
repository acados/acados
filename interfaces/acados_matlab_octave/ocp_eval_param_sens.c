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
#include "acados_c/ocp_nlp_interface.h"
// mex
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    long long *ptr;

    /* RHS */

    // C_ocp

    // solver
    ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "solver" ) );
    ocp_nlp_solver *solver = (ocp_nlp_solver *) ptr[0];
    // sens_out
    ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "sens_out" ) );
    ocp_nlp_out *sens_out = (ocp_nlp_out *) ptr[0];

    // field
    char *field = mxArrayToString( prhs[1] );

    // stage
    int stage = mxGetScalar( prhs[2] );

    // index
    int index = mxGetScalar( prhs[3] );

    /* solver */
    ocp_nlp_eval_param_sens(solver, field, stage, index, sens_out);


    return;

}



