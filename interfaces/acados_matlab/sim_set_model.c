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
	char *method = mxArrayToString( mxGetField( prhs[0], 0, "method" ) );
//	mexPrintf("\n%s\n", method);


	// C_sim

	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "config" ) );
	sim_config *config = (sim_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dims" ) );
	void *dims = (void *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "in" ) );
	sim_in *in = (sim_in *) ptr[0];


	// C_sim_ext_fun

	
	// TODO check for empty struct member


	/* set in model */

	if(!strcmp(method, "erk"))
		{
//		mexPrintf("\n%s\n", method);

		// expl_ode_fun
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_expl_ode_fun" ) );
		external_function_casadi *expl_ode_fun = (external_function_casadi *) ptr[0];
		// expl_vde_for
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_expl_vde_for" ) );
		external_function_casadi *expl_vde_for = (external_function_casadi *) ptr[0];
		// expl_vde_adj
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_expl_vde_adj" ) );
		external_function_casadi *expl_vde_adj = (external_function_casadi *) ptr[0];

		sim_in_set(config, dims, in, "expl_ode_fun", expl_ode_fun);
		sim_in_set(config, dims, in, "expl_vde_for", expl_vde_for);
		}
	else if(!strcmp(method, "irk"))
		{
//		mexPrintf("\n%s\n", method);

		// impl_ode_fun
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_impl_ode_fun" ) );
		external_function_casadi *impl_ode_fun = (external_function_casadi *) ptr[0];
		// impl_ode_fun_jac_x_xdot
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_impl_ode_fun_jac_x_xdot" ) );
		external_function_casadi *impl_ode_fun_jac_x_xdot = (external_function_casadi *) ptr[0];
		// impl_ode_jac_x_xdot_u
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_impl_ode_jac_x_xdot_u" ) );
		external_function_casadi *impl_ode_jac_x_xdot_u = (external_function_casadi *) ptr[0];

		sim_in_set(config, dims, in, "impl_ode_fun", impl_ode_fun);
		sim_in_set(config, dims, in, "impl_ode_fun_jac_x_xdot", impl_ode_fun_jac_x_xdot);
		sim_in_set(config, dims, in, "impl_ode_jac_x_xdot_u", impl_ode_jac_x_xdot_u);
		}
	else
		{
		mexPrintf("\nsim_set_model: method not supported %s\n", method);
		return;
		}



	/* return */

	return;

	}

