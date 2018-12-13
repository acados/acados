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

#if 0
	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "config" ) );
	ocp_nlp_config *config = (sim_config *) ptr[0];
//	mexPrintf("\n%lld %p\n", config_mat[0], config_mat[0]);
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "dims" ) );
	void *dims = (void *) ptr[0];
	// opts
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "opts" ) );
	sim_opts *opts = (sim_opts *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "in" ) );
	sim_in *in = (sim_in *) ptr[0];
	// out
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "out" ) );
	sim_out *out = (sim_out *) ptr[0];
	// solver
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "solver" ) );
	sim_solver *solver = (sim_solver *) ptr[0];



	/* free memory */

	ocp_nlp_config_destroy(config);
	ocp_nlp_dims_destroy(dims);
	ocp_nlp_opts_destroy(opts);
	ocp_nlp_in_destroy(in);
	ocp_nlp_out_destroy(out);
	ocp_nlp_solver_destroy(solver);
#endif



	/* return */

	return;

	}

