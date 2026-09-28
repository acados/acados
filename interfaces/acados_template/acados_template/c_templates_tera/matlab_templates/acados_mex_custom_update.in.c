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
#include "acados_solver_{{ name }}.h"

// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // C_ocp
    long long *ptr;
    const mxArray *C_ocp = prhs[0];

    // capsule
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "capsule" ) );
    {{ name }}_solver_capsule *capsule = ({{ name }}_solver_capsule *) ptr[0];

    // data
    double *data = mxGetPr(prhs[1]);
    int data_len = (int) mxGetNumberOfElements( prhs[1] );

    // solve
    int status = {{ name }}_acados_custom_update(capsule, data, data_len);

    plhs[0] = mxCreateNumericMatrix(1, 11, mxDOUBLE_CLASS, mxREAL);
    double *x = mxGetPr( plhs[0] );
    x[0] = (double) status;

    return;
}
