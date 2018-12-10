// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados/sim/sim_common.h"
#include "acados_c/sim_interface.h"
#include "acados/utils/external_function_generic.h"
#include "acados_c/external_function_interface.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin sim_set_model\n");

	long long *ptr;


	/* RHS */

	// opts_struct

	int num_stages = mxGetScalar( mxGetField( prhs[0], 0, "num_stages" ) );
	int num_steps = mxGetScalar( mxGetField( prhs[0], 0, "num_steps" ) );
	bool sens_forw = mxGetScalar( mxGetField( prhs[0], 0, "sens_forw" ) );
//	mexPrintf("\n%d\n", sens_forw);
	char *scheme = mxArrayToString( mxGetField( prhs[0], 0, "scheme" ) );
//	mexPrintf("\n%s\n", scheme);


	// C_sim

	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "config" ) );
	sim_solver_config *config = (sim_solver_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dims" ) );
	void *dims = (void *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "in" ) );
	sim_in *in = (sim_in *) ptr[0];


	// C_sim_ext_fun

	// expl_ode_fun
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "expl_ode_fun" ) );
	external_function_casadi *expl_ode_fun = (external_function_casadi *) ptr[0];
	// expl_vde_for
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "expl_vde_for" ) );
	external_function_casadi *expl_vde_for = (external_function_casadi *) ptr[0];



	/* set in model */

	if(!strcmp(scheme, "erk"))
		{
//		mexPrintf("\n%s\n", scheme);
		sim_in_set(config, dims, in, "expl_ode_fun", expl_ode_fun);
		sim_in_set(config, dims, in, "expl_vde_for", expl_vde_for);
		}
	else
		{
		mexPrintf("\nsim_set_model: scheme not supported %s\n", scheme);
		return;
		}



	/* return */

	return;

	}

