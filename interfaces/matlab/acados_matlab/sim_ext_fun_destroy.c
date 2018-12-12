// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados/utils/external_function_generic.h"
#include "acados_c/external_function_interface.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin sim_ext_fun_destroy\n");

	long long *ptr;



	/* RHS */

	// opts_struct

	bool sens_forw = mxGetScalar( mxGetField( prhs[0], 0, "sens_forw" ) );
//	mexPrintf("\n%d\n", sens_forw);
	char *scheme = mxArrayToString( mxGetField( prhs[0], 0, "scheme" ) );
//	mexPrintf("\n%s\n", scheme);


	// TODO check for empty struct member


	/* free memory */

	if(!strcmp(scheme, "erk"))
		{
		// C_sim_ext_fun

		// expl_ode_fun
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "expl_ode_fun" ) );
		external_function_casadi *expl_ode_fun = (external_function_casadi *) ptr[0];
		// expl_vde_for
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "expl_vde_for" ) );
		external_function_casadi *expl_vde_for = (external_function_casadi *) ptr[0];
		// expl_vde_adj
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "expl_vde_adj" ) );
		external_function_casadi *expl_vde_adj = (external_function_casadi *) ptr[0];

		// free external functions
		external_function_casadi_free(expl_ode_fun);
		external_function_casadi_free(expl_vde_for);
		free(expl_ode_fun);
		free(expl_vde_for);
		}
	else if(!strcmp(scheme, "irk"))
		{
		// C_sim_ext_fun

		// impl_ode_fun
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "impl_ode_fun" ) );
		external_function_casadi *impl_ode_fun = (external_function_casadi *) ptr[0];
		// impl_ode_fun_jac_x_xdot
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "impl_ode_fun_jac_x_xdot" ) );
		external_function_casadi *impl_ode_fun_jac_x_xdot = (external_function_casadi *) ptr[0];
		// impl_ode_jac_x_xdot_u
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "impl_ode_jac_x_xdot_u" ) );

		// free external functions
		external_function_casadi *impl_ode_jac_x_xdot_u = (external_function_casadi *) ptr[0];
		external_function_casadi_free(impl_ode_fun);
		external_function_casadi_free(impl_ode_fun_jac_x_xdot);
		external_function_casadi_free(impl_ode_jac_x_xdot_u);
		free(impl_ode_fun);
		free(impl_ode_fun_jac_x_xdot);
		free(impl_ode_jac_x_xdot_u);
		}
	else
		{
		mexPrintf("\nscheme not supported %s\n", scheme);
		return;
		}




	/* return */

	return;

	}

