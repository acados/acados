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

//	mexPrintf("\nin sim_destroy\n");

	long long *ptr;

//	void *config = mxGetPr( mxGetField( prhs[0], 0, "config" ) );
//	long long *config_mat = (long long *) mxGetData( mxGetField( prhs[0], 0, "config" ) );
//	long long config = (long long) mxGetScalar( mxGetField( prhs[0], 0, "config" ) );



	/* RHS */

	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "config" ) );
	sim_config *config = (sim_config *) ptr[0];
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

	sim_config_destroy(config);
	sim_dims_destroy(dims);
	sim_opts_destroy(opts);
	sim_in_destroy(in);
	sim_out_destroy(out);
	sim_solver_destroy(solver);



	/* return */

	return;

	}
