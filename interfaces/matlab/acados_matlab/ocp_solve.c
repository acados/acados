// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
//#include "acados/sim/sim_common.h"
#include "acados_c/ocp_nlp_interface.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin ocp_solve\n");

	long long *ptr;

	/* RHS */

	// C_ocp

	// solver
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "solver" ) );
	ocp_nlp_solver *solver = (ocp_nlp_solver *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "in" ) );
	ocp_nlp_in *in = (ocp_nlp_in *) ptr[0];
	// out
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "out" ) );
	ocp_nlp_out *out = (ocp_nlp_out *) ptr[0];



	/* solver */
	int acados_return = ocp_nlp_solve(solver, in, out);



	/* return */

	return;

	}



