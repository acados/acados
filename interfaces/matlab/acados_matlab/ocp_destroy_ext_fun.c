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

//	mexPrintf("\nin ocp_ext_fun_destroy\n");

	long long *ptr;



	/* RHS */

	// C_ocp_ext_fun

	external_function_casadi *ext_fun_ptr;

	if (mxGetField( prhs[1], 0, "expl_ode_fun" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "expl_ode_fun" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		external_function_casadi_free(ext_fun_ptr);
		free(ext_fun_ptr);
		}
	if (mxGetField( prhs[1], 0, "expl_vde_for" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "expl_vde_for" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		external_function_casadi_free(ext_fun_ptr);
		free(ext_fun_ptr);
		}
	if (mxGetField( prhs[1], 0, "expl_vde_adj" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "expl_vde_adj" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		external_function_casadi_free(ext_fun_ptr);
		free(ext_fun_ptr);
		}
	if (mxGetField( prhs[1], 0, "impl_ode_fun" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "impl_ode_fun" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		external_function_casadi_free(ext_fun_ptr);
		free(ext_fun_ptr);
		}
	if (mxGetField( prhs[1], 0, "impl_ode_fun_jac_x_xdot" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "impl_ode_fun_jac_x_xdot" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		external_function_casadi_free(ext_fun_ptr);
		free(ext_fun_ptr);
		}
	if (mxGetField( prhs[1], 0, "impl_ode_jac_x_xdot_u" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "impl_ode_jac_x_xdot_u" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		external_function_casadi_free(ext_fun_ptr);
		free(ext_fun_ptr);
		}
	if (mxGetField( prhs[1], 0, "h_fun_jac_ut_xt" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "h_fun_jac_ut_xt" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		external_function_casadi_free(ext_fun_ptr);
		free(ext_fun_ptr);
		}



	/* return */

	return;

	}


