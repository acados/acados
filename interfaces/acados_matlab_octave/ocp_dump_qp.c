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

extern void ocp_nlp_dump_last_qp_to_json(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_solver *solver, const char *filename);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    long long *ptr;
    char filename[256];

    // C_ocp
    const mxArray *C_ocp = prhs[0];
    // config
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "config" ) );
    ocp_nlp_config *config = (ocp_nlp_config *) ptr[0];
    // dims
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "dims" ) );
    ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];
    // solver
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "solver" ) );
    ocp_nlp_solver *solver = (ocp_nlp_solver *) ptr[0];

    mxGetString(prhs[1], filename, sizeof(filename));

    ocp_nlp_dump_last_qp_to_json(config, dims, solver, filename);

    return;
}