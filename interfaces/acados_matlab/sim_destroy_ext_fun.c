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

	// model

	// TODO bool instead !!!
	char *param_f = mxArrayToString( mxGetField( prhs[0], 0, "dyn_param_f" ) );

	// C_sim_ext_fun

	external_function_casadi *ext_fun_ptr;
	external_function_param_casadi *ext_fun_param_ptr;

	if(mxGetField( prhs[1], 0, "dyn_expl_ode_fun" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dyn_expl_ode_fun" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			external_function_param_casadi_free(ext_fun_param_ptr);
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			external_function_casadi_free(ext_fun_ptr);
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[1], 0, "dyn_expl_vde_for" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dyn_expl_vde_for" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			external_function_param_casadi_free(ext_fun_param_ptr);
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			external_function_casadi_free(ext_fun_ptr);
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[1], 0, "dyn_expl_vde_adj" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dyn_expl_vde_adj" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			external_function_param_casadi_free(ext_fun_param_ptr);
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			external_function_casadi_free(ext_fun_ptr);
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[1], 0, "dyn_impl_ode_fun" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dyn_impl_ode_fun" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			external_function_param_casadi_free(ext_fun_param_ptr);
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			external_function_casadi_free(ext_fun_ptr);
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[1], 0, "dyn_impl_ode_fun_jac_x_xdot" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dyn_impl_ode_fun_jac_x_xdot" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			external_function_param_casadi_free(ext_fun_param_ptr);
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			external_function_casadi_free(ext_fun_ptr);
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[1], 0, "dyn_impl_ode_jac_x_xdot_u" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dyn_impl_ode_jac_x_xdot_u" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			external_function_param_casadi_free(ext_fun_param_ptr);
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			external_function_casadi_free(ext_fun_ptr);
			free(ext_fun_ptr);
			}
		}



	/* return */

	return;

	}

