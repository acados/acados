// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados/sim/sim_common.h"
#include "acados_c/sim_interface.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin sim_set\n");

	long long *ptr;

	/* RHS */

	// C_sim

	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "config" ) );
	sim_config *config = (sim_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "dims" ) );
	void *dims = (void *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "in" ) );
	sim_in *in = (sim_in *) ptr[0];

	// field
	char *field = mxArrayToString( prhs[1] );
//	mexPrintf("\n%s\n", field);


	// value
	if(!strcmp(field, "T"))
		{
		double *T = mxGetPr( prhs[2] );
		sim_in_set(config, dims, in, "T", T);
		}
	else if(!strcmp(field, "x"))
		{
		double *x = mxGetPr( prhs[2] );
		sim_in_set(config, dims, in, "x", x);
		}
	else if(!strcmp(field, "u"))
		{
		double *u = mxGetPr( prhs[2] );
		sim_in_set(config, dims, in, "u", u);
		}
	else
		{
		mexPrintf("\nsim_set: field not supported: %s\n", field);
		return;
		}



	/* return */

	return;

	}


