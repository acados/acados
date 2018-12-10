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



#if 0
	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "config" ) );
	sim_solver_config *config = (sim_solver_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "dims" ) );
	void *dims = (void *) ptr[0];
	int nx;
	config->get_nx(dims, &nx);
	mexPrintf("\nnx %d\n", nx);
	d_print_mat(1, nx, out->xn, 1);
#endif



	/* return */

	return;

	}


