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
// example specific
#include "acados_sim_solver_{{ name }}.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* RHS */
    const mxArray *C_sim = prhs[0];
    long long * ptr;

    // capsule
    ptr = (long long *) mxGetData( mxGetField( C_sim, 0, "capsule" ) );
    {{ name }}_sim_solver_capsule *capsule = ({{ name }}_sim_solver_capsule *) ptr[0];


    /* free memory */
    int status = 0;

    status = {{ name }}_acados_sim_free(capsule);
    if (status)
    {
        mexPrintf("{{ name }}_acados_sim_free() returned status %d.\n", status);
    }

    status = {{ name }}_acados_sim_solver_free_capsule(capsule);
    if (status)
    {
        mexPrintf("{{ name }}_acados_sim_solver_free_capsule() returned status %d.\n", status);
    }

    return;

}
