/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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
#include "acados_c/ocp_nlp_interface.h"
#include "acados/dense_qp/dense_qp_common.h"
// mex
#include "mex.h"
#include "mex_macros.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    long long *ptr;
    char fun_name[50] = "ocp_get";
    char buffer [200]; // for error messages

    int ii, jj;

    /* RHS */

    // C_ocp
    const mxArray *C_ocp = prhs[0];
    // config
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "config" ) );
    ocp_nlp_config *config = (ocp_nlp_config *) ptr[0];
    // dims
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "dims" ) );
    ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];
    // out
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "out" ) );
    ocp_nlp_out *out = (ocp_nlp_out *) ptr[0];
    // solver
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "solver" ) );
    ocp_nlp_solver *solver = (ocp_nlp_solver *) ptr[0];
    // sens_out
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "sens_out" ) );
    ocp_nlp_out *sens_out = (ocp_nlp_out *) ptr[0];
    // plan
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "plan" ) );
    ocp_nlp_plan *plan = (ocp_nlp_plan *) ptr[0];

    // field
    char *field = mxArrayToString( prhs[1] );
//    mexPrintf("\nin ocp_get: field%s\n", field);

    int N = dims->N;
    int nu = dims->nu[0];
    int nx = dims->nx[0];
    int nz = dims->nz[0];

    if (!strcmp(field, "x"))
    {
        if (nrhs==2)
        {
            plhs[0] = mxCreateNumericMatrix(nx, N+1, mxDOUBLE_CLASS, mxREAL);
            double *x = mxGetPr( plhs[0] );
            for (ii=0; ii<=N; ii++)
            {
                ocp_nlp_out_get(config, dims, out, ii, "x", x+ii*nx);
            }
        }
        else if (nrhs==3)
        {
            plhs[0] = mxCreateNumericMatrix(nx, 1, mxDOUBLE_CLASS, mxREAL);
            double *x = mxGetPr( plhs[0] );
            int stage = mxGetScalar( prhs[2] );
            ocp_nlp_out_get(config, dims, out, stage, "x", x);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "u"))
    {
        if (nrhs==2)
        {
            plhs[0] = mxCreateNumericMatrix(nu, N, mxDOUBLE_CLASS, mxREAL);
            double *u = mxGetPr( plhs[0] );
            for (ii=0; ii<N; ii++)
            {
                ocp_nlp_out_get(config, dims, out, ii, "u", u+ii*nu);
            }
        }
        else if (nrhs==3)
        {
            plhs[0] = mxCreateNumericMatrix(nu, 1, mxDOUBLE_CLASS, mxREAL);
            double *u = mxGetPr( plhs[0] );
            int stage = mxGetScalar( prhs[2] );
            ocp_nlp_out_get(config, dims, out, stage, "u", u);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "z"))
    {
        if (nrhs==2)
        {
            plhs[0] = mxCreateNumericMatrix(nz, N, mxDOUBLE_CLASS, mxREAL);
            double *z = mxGetPr( plhs[0] );
            for (ii=0; ii<N; ii++)
            {
                ocp_nlp_out_get(config, dims, out, ii, "z", z+ii*nz);
            }
        }
        else if (nrhs==3)
        {
            plhs[0] = mxCreateNumericMatrix(nz, 1, mxDOUBLE_CLASS, mxREAL);
            double *z = mxGetPr( plhs[0] );
            int stage = mxGetScalar( prhs[2] );
            ocp_nlp_out_get(config, dims, out, stage, "z", z);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "pi"))
    {
        if (nrhs==2)
        {
            plhs[0] = mxCreateNumericMatrix(nx, N, mxDOUBLE_CLASS, mxREAL);
            double *pi = mxGetPr( plhs[0] );
            for (ii=0; ii<N; ii++)
            {
                ocp_nlp_out_get(config, dims, out, ii, "pi", pi+ii*nx);
            }
        }
        else if (nrhs==3)
        {
            plhs[0] = mxCreateNumericMatrix(nx, 1, mxDOUBLE_CLASS, mxREAL);
            double *pi = mxGetPr( plhs[0] );
            int stage = mxGetScalar( prhs[2] );
            ocp_nlp_out_get(config, dims, out, stage, "pi", pi);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "sens_x"))
    {
        if (nrhs==2)
        {
            plhs[0] = mxCreateNumericMatrix(nx, N+1, mxDOUBLE_CLASS, mxREAL);
            double *x = mxGetPr( plhs[0] );
            for (ii=0; ii<=N; ii++)
            {
                ocp_nlp_out_get(config, dims, sens_out, ii, "x", x+ii*nx);
            }
        }
        else if (nrhs==3)
        {
            plhs[0] = mxCreateNumericMatrix(nx, 1, mxDOUBLE_CLASS, mxREAL);
            double *x = mxGetPr( plhs[0] );
            int stage = mxGetScalar( prhs[2] );
            ocp_nlp_out_get(config, dims, sens_out, stage, "x", x);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "sens_u"))
    {
        if (nrhs==2)
        {
            plhs[0] = mxCreateNumericMatrix(nu, N, mxDOUBLE_CLASS, mxREAL);
            double *u = mxGetPr( plhs[0] );
            for (ii=0; ii<N; ii++)
            {
                ocp_nlp_out_get(config, dims, sens_out, ii, "u", u+ii*nu);
            }
        }
        else if (nrhs==3)
        {
            plhs[0] = mxCreateNumericMatrix(nu, 1, mxDOUBLE_CLASS, mxREAL);
            double *u = mxGetPr( plhs[0] );
            int stage = mxGetScalar( prhs[2] );
            ocp_nlp_out_get(config, dims, sens_out, stage, "u", u);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "sens_pi"))
    {
        if (nrhs==2)
        {
            plhs[0] = mxCreateNumericMatrix(nx, N, mxDOUBLE_CLASS, mxREAL);
            double *pi = mxGetPr( plhs[0] );
            for (ii=0; ii<N; ii++)
            {
                ocp_nlp_out_get(config, dims, sens_out, ii, "pi", pi+ii*nx);
            }
        }
        else if (nrhs==3)
        {
            plhs[0] = mxCreateNumericMatrix(nx, 1, mxDOUBLE_CLASS, mxREAL);
            double *pi = mxGetPr( plhs[0] );
            int stage = mxGetScalar( prhs[2] );
            ocp_nlp_out_get(config, dims, sens_out, stage, "pi", pi);
        }
        else
        {
            mexPrintf("\nocp_get: wrong nrhs: %d\n", nrhs);
            return;
        }
    }
    else if (!strcmp(field, "status"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        int status;
        ocp_nlp_get(config, solver, "status", &status);
        *mat_ptr = (double) status;
    }
    else if (!strcmp(field, "sqp_iter"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        int sqp_iter;
        ocp_nlp_get(config, solver, "sqp_iter", &sqp_iter);
        *mat_ptr = (double) sqp_iter;
    }
    else if (!strcmp(field, "time_tot"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        ocp_nlp_get(config, solver, "time_tot", mat_ptr);
    }
    else if (!strcmp(field, "time_lin"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        ocp_nlp_get(config, solver, "time_lin", mat_ptr);
    }
    else if (!strcmp(field, "time_reg"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        ocp_nlp_get(config, solver, "time_reg", mat_ptr);
    }
    else if (!strcmp(field, "time_qp_sol"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        ocp_nlp_get(config, solver, "time_qp_sol", mat_ptr);
    }
    else if (!strcmp(field, "stat"))
    {
        int sqp_iter;
        int stat_m, stat_n;
        double *stat;
        ocp_nlp_get(config, solver, "sqp_iter", &sqp_iter);
        ocp_nlp_get(config, solver, "stat_m", &stat_m);
        ocp_nlp_get(config, solver, "stat_n", &stat_n);
        ocp_nlp_get(config, solver, "stat", &stat);
        int min_size = stat_m<sqp_iter+1 ? stat_m : sqp_iter+1;
        plhs[0] = mxCreateNumericMatrix(min_size, stat_n+1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        for (ii=0; ii<min_size; ii++)
        {
            mat_ptr[ii+0] = ii;
            for (jj=0; jj<stat_n; jj++)
                mat_ptr[ii+(jj+1)*min_size] = stat[jj+ii*stat_n];
        }
    }
    else if (!strcmp(field, "qp_solver_hessian"))
    {
		void *qp_in_;
        ocp_nlp_get(config, solver, "qp_xcond_in", &qp_in_);
		int solver_type = 0;
		if (plan->ocp_qp_solver_plan.qp_solver==FULL_CONDENSING_HPIPM)
			solver_type=1;
#if defined(ACADOS_WITH_QPOASES)
		if (plan->ocp_qp_solver_plan.qp_solver==FULL_CONDENSING_QPOASES)
			solver_type=1;
#endif
		// dense solver
		if(solver_type==1)
		{
			dense_qp_in *qp_in = qp_in_;
			int nv = qp_in->dim->nv;
			plhs[0] = mxCreateNumericMatrix(nv, nv, mxDOUBLE_CLASS, mxREAL);
			double *mat_ptr = mxGetPr( plhs[0] );
			blasfeo_unpack_dmat(nv, nv, qp_in->Hv, 0, 0, mat_ptr, nv);
		}
		// TODO ocp qp hessian
		else
		{
			plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
		}
	}
    else
    {
        MEX_FIELD_NOT_SUPPORTED_SUGGEST(fun_name, field,
             "x, u, z, pi, sens_x, sens_u, sens_pi, status, sqp_iter, time_tot, time_lin, time_reg, time_qp_sol, stat");
    }

    return;

}

