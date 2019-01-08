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

//	mexPrintf("\nin sim_solve\n");

	long long *ptr;

	/* RHS */

	// C_sim

	// solver
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "solver" ) );
	sim_solver *solver = (sim_solver *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "in" ) );
	sim_in *in = (sim_in *) ptr[0];
	// out
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "out" ) );
	sim_out *out = (sim_out *) ptr[0];



	/* solver */
	int acados_return = sim_solve(solver, in, out);



	/* return */

	return;

	}


