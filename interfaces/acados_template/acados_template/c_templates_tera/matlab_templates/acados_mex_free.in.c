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
    int status = 0;
    long long *ptr;

    // mexPrintf("\nin mex_acados_free\n");
    const mxArray *C_ocp = prhs[0];
    // capsule
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "capsule" ) );
    {{ name }}_solver_capsule *capsule = ({{ name }}_solver_capsule *) ptr[0];

    status = {{ name }}_acados_free(capsule);
    if (status)
    {
        mexPrintf("{{ name }}_acados_free() returned status %d.\n", status);
    }

    status = {{ name }}_acados_free_capsule(capsule);
    if (status)
    {
        mexPrintf("{{ name }}_acados_free_capsule() returned status %d.\n", status);
    }

    return;
}

