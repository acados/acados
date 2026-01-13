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
#include "acados_solver_{{ name }}.h"

// mex
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Usage:
    //   P_all = acados_mex_get_zoRO_Pk_{{ name }}(C_ocp)           -> returns 1x(N+1) cell
    //   P_k   = acados_mex_get_zoRO_Pk_{{ name }}(C_ocp, stage)    -> returns [nx x nx] matrix

    if (nrhs < 1 || nrhs > 2)
    {
        mexErrMsgTxt("acados_mex_get_zoRO_Pk:  Expected 1 or 2 arguments:  (C_ocp) or (C_ocp, stage)");
    }

    // C_ocp
    long long *ptr;
    const mxArray *C_ocp = prhs[0];

    // capsule
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "capsule" ) );
    {{ name }}_solver_capsule *capsule = ({{ name }}_solver_capsule *) ptr[0];

    // Get dimensions
    ocp_nlp_dims *dims = {{ name }}_acados_get_nlp_dims(capsule);
    int N = dims->N;
    int nx = dims->nx[0];
    int total_size = (N + 1) * nx * nx;

    // Allocate buffer and get all P matrices
    double *P_buffer = (double *) mxMalloc(total_size * sizeof(double));
    int status = {{ name }}_acados_get_zoRO_Pk_matrices(capsule, P_buffer, total_size);

    if (status != 0)
    {
        mxFree(P_buffer);
        mexErrMsgTxt("acados_mex_get_zoRO_Pk: Failed to get zoRO P matrices");
    }

    if (nrhs == 1)
    {
        // Return cell array of all P matrices
        plhs[0] = mxCreateCellMatrix(1, N + 1);
        for (int k = 0; k <= N; k++)
        {
            mxArray *Pk = mxCreateDoubleMatrix(nx, nx, mxREAL);
            memcpy(mxGetPr(Pk), &P_buffer[k * nx * nx], nx * nx * sizeof(double));
            mxSetCell(plhs[0], k, Pk);
        }
    }
    else
    {
        // Return single P matrix for requested stage
        int stage = (int) mxGetScalar(prhs[1]);
        if (stage < 0 || stage > N)
        {
            mxFree(P_buffer);
            mexErrMsgTxt("acados_mex_get_zoRO_Pk:  stage out of bounds");
        }
        plhs[0] = mxCreateDoubleMatrix(nx, nx, mxREAL);
        memcpy(mxGetPr(plhs[0]), &P_buffer[stage * nx * nx], nx * nx * sizeof(double));
    }

    mxFree(P_buffer);
    return;
}