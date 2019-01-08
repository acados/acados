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

//	mexPrintf("\nin sim_destroy\n");

	long long *ptr;

//	void *config = mxGetPr( mxGetField( prhs[0], 0, "config" ) );
//	long long *config_mat = (long long *) mxGetData( mxGetField( prhs[0], 0, "config" ) );
//	long long config = (long long) mxGetScalar( mxGetField( prhs[0], 0, "config" ) );



	/* RHS */

	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "config" ) );
	ocp_nlp_config *config = (ocp_nlp_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "dims" ) );
	ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];
	// opts
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "opts" ) );
	void *opts = (void *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "in" ) );
	ocp_nlp_in *in = (ocp_nlp_in *) ptr[0];
	// out
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "out" ) );
	ocp_nlp_out *out = (ocp_nlp_out *) ptr[0];
	// solver
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "solver" ) );
	ocp_nlp_solver *solver = (ocp_nlp_solver *) ptr[0];



	/* free memory */

	ocp_nlp_config_destroy(config);
	ocp_nlp_dims_destroy(dims);
	ocp_nlp_opts_destroy(opts);
	ocp_nlp_in_destroy(in);
	ocp_nlp_out_destroy(out);
	ocp_nlp_solver_destroy(solver);



	/* return */

	return;

	}

