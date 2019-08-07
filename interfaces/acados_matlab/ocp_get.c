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

//	mexPrintf("\nin ocp_get\n");

	long long *ptr;

	int ii, jj;

	/* RHS */

	// C_ocp

	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "config" ) );
	ocp_nlp_config *config = (ocp_nlp_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "dims" ) );
	ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];
	// out
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "out" ) );
	ocp_nlp_out *out = (ocp_nlp_out *) ptr[0];
	// solver
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "solver" ) );
	ocp_nlp_solver *solver = (ocp_nlp_solver *) ptr[0];

	// field
	char *field = mxArrayToString( prhs[1] );
//	mexPrintf("\n%s\n", field);



	int N = dims->N;
	int nu = dims->nu[0];
	int nx = dims->nx[0];



	if(!strcmp(field, "x"))
		{
		if(nrhs==2)
			{
			plhs[0] = mxCreateNumericMatrix(nx, N+1, mxDOUBLE_CLASS, mxREAL);
			double *x = mxGetPr( plhs[0] );
			for(ii=0; ii<=N; ii++)
				{
				ocp_nlp_out_get(config, dims, out, ii, "x", x+ii*nx);
				}
			}
		else if(nrhs==3)
			{
			plhs[0] = mxCreateNumericMatrix(nx, 1, mxDOUBLE_CLASS, mxREAL);
			double *x = mxGetPr( plhs[0] );
			int stage = mxGetScalar( prhs[2] );
			ocp_nlp_out_get(config, dims, out, stage, "x", x);
			}
		else
			{
			mexPrintf("\nocp_get: wrong nrhs: %d\n", nrhs);
			return;
			}
		}
	else if(!strcmp(field, "u"))
		{
		if(nrhs==2)
			{
			plhs[0] = mxCreateNumericMatrix(nu, N, mxDOUBLE_CLASS, mxREAL);
			double *u = mxGetPr( plhs[0] );
			for(ii=0; ii<N; ii++)
				{
				ocp_nlp_out_get(config, dims, out, ii, "u", u+ii*nu);
				}
			}
		else if(nrhs==3)
			{
			plhs[0] = mxCreateNumericMatrix(nu, 1, mxDOUBLE_CLASS, mxREAL);
			double *u = mxGetPr( plhs[0] );
			int stage = mxGetScalar( prhs[2] );
			ocp_nlp_out_get(config, dims, out, stage, "u", u);
			}
		else
			{
			mexPrintf("\nocp_get: wrong nrhs: %d\n", nrhs);
			return;
			}
		}
	else if(!strcmp(field, "pi"))
		{
		if(nrhs==2)
			{
			plhs[0] = mxCreateNumericMatrix(nx, N, mxDOUBLE_CLASS, mxREAL);
			double *pi = mxGetPr( plhs[0] );
			for(ii=0; ii<N; ii++)
				{
				ocp_nlp_out_get(config, dims, out, ii, "pi", pi+ii*nx);
				}
			}
		else if(nrhs==3)
			{
			plhs[0] = mxCreateNumericMatrix(nx, 1, mxDOUBLE_CLASS, mxREAL);
			double *pi = mxGetPr( plhs[0] );
			int stage = mxGetScalar( prhs[2] );
			ocp_nlp_out_get(config, dims, out, stage, "pi", pi);
			}
		else
			{
			mexPrintf("\nocp_get: wrong nrhs: %d\n", nrhs);
			return;
			}
		}
	else if(!strcmp(field, "status"))
		{
		plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		double *mat_ptr = mxGetPr( plhs[0] );
		int status;
		ocp_nlp_get(config, solver, "status", &status);
		*mat_ptr = (double) status;
		}
	else if(!strcmp(field, "sqp_iter"))
		{
		plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		double *mat_ptr = mxGetPr( plhs[0] );
		int sqp_iter;
		ocp_nlp_get(config, solver, "sqp_iter", &sqp_iter);
		*mat_ptr = (double) sqp_iter;
		}
	else if(!strcmp(field, "time_tot"))
		{
		plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		double *mat_ptr = mxGetPr( plhs[0] );
		ocp_nlp_get(config, solver, "time_tot", mat_ptr);
		}
	else if(!strcmp(field, "time_lin"))
		{
		plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		double *mat_ptr = mxGetPr( plhs[0] );
		ocp_nlp_get(config, solver, "time_lin", mat_ptr);
		}
	else if(!strcmp(field, "time_reg"))
		{
		plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		double *mat_ptr = mxGetPr( plhs[0] );
		ocp_nlp_get(config, solver, "time_reg", mat_ptr);
		}
	else if(!strcmp(field, "time_qp_sol"))
		{
		plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		double *mat_ptr = mxGetPr( plhs[0] );
		ocp_nlp_get(config, solver, "time_qp_sol", mat_ptr);
		}
	else if(!strcmp(field, "stat"))
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
		for(ii=0; ii<min_size; ii++)
			{
			mat_ptr[ii+0] = ii;
			for(jj=0; jj<stat_n; jj++)
				mat_ptr[ii+(jj+1)*min_size] = stat[jj+ii*stat_n];
			}
		}
	else
		{
		mexPrintf("\nocp_get: field not supported: %s\n", field);
		return;
		}



	/* return */

	return;

	}



