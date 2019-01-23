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
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "in" ) );
	ocp_nlp_in *in = (ocp_nlp_in *) ptr[0];
	// out
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "out" ) );
	ocp_nlp_out *out = (ocp_nlp_out *) ptr[0];

	// field
	char *field = mxArrayToString( prhs[1] );
//	mexPrintf("\n%s\n", field);



	int N = dims->N;
	int nu = dims->nu[0];
	int nx = dims->nx[0];



	// TODO implement with LHS !!!!!
	// value
	if (!strcmp(field, "x0"))
		{
		double *x0 = mxGetPr( prhs[2] );
		ocp_nlp_constraints_model_set(config, dims, in, 0, "lbx", x0);
		ocp_nlp_constraints_model_set(config, dims, in, 0, "ubx", x0);
		}
	else if (!strcmp(field, "x_init"))
		{
		double *x_init = mxGetPr( prhs[2] );
		for (ii=0; ii<=N; ii++)
			{
			ocp_nlp_out_set(config, dims, out, ii, "x", x_init+ii*nx);
			}
		}
	else if (!strcmp(field, "u_init"))
		{
		double *u_init = mxGetPr( prhs[2] );
		for (ii=0; ii<N; ii++)
			{
			ocp_nlp_out_set(config, dims, out, ii, "u", u_init+ii*nu);
			}
		}
	// TODO make setters for ALL numerical data
	else
		{
		mexPrintf("\nocp_set: field not supported: %s\n", field);
		return;
		}



	/* return */

	return;

	}




