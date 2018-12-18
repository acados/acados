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

	int ii;

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

	// field
	char *field = mxArrayToString( prhs[1] );
//	mexPrintf("\n%s\n", field);



	int N = dims->N;
	int nu = dims->nu[0];
	int nx = dims->nx[0];



	if(!strcmp(field, "x"))
		{
		plhs[0] = mxCreateNumericMatrix(nx, N+1, mxDOUBLE_CLASS, mxREAL);
		double *x = mxGetPr( plhs[0] );
		for(ii=0; ii<=N; ii++)
			ocp_nlp_out_get(config, dims, out, ii, "x", x+ii*nx);
		}
	else if(!strcmp(field, "u"))
		{
		plhs[0] = mxCreateNumericMatrix(nu, N, mxDOUBLE_CLASS, mxREAL);
		double *u = mxGetPr( plhs[0] );
		for(ii=0; ii<N; ii++)
			ocp_nlp_out_get(config, dims, out, ii, "u", u+ii*nu);
		}
	else
		{
		mexPrintf("\nocp_get: field not supported: %s\n", field);
		return;
		}



	/* return */

	return;

	}



