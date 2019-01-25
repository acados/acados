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

//	mexPrintf("\nin sim_set\n");

	long long *ptr;

	/* RHS */

	// model
	// TODO bool instead !!!
	char *param = mxArrayToString( mxGetField( prhs[0], 0, "param" ) );

	// opts
	char *method = mxArrayToString( mxGetField( prhs[1], 0, "method" ) );

	// C_sim

	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "config" ) );
	sim_config *config = (sim_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dims" ) );
	void *dims = (void *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "in" ) );
	sim_in *in = (sim_in *) ptr[0];

	// field
	char *field = mxArrayToString( prhs[4] );
//	mexPrintf("\n%s\n", field);


	// value
	if(!strcmp(field, "T"))
		{
		double *T = mxGetPr( prhs[5] );
		sim_in_set(config, dims, in, "T", T);
		}
	else if(!strcmp(field, "x"))
		{
		double *x = mxGetPr( prhs[5] );
		sim_in_set(config, dims, in, "x", x);
		}
	else if(!strcmp(field, "u"))
		{
		double *u = mxGetPr( prhs[5] );
		sim_in_set(config, dims, in, "u", u);
		}
	else if(!strcmp(field, "p"))
		{
		if(!strcmp(param, "true")) // TODO bool
			{
			double *p = mxGetPr( prhs[5] );
			// TODO check for param in opts or the like !!!!!!!!!!!!!!!!!!!!!!!
			if(!strcmp(method, "erk"))
				{
				// expl_ode_fun
				ptr = (long long *) mxGetData( mxGetField( prhs[3], 0, "expl_ode_fun" ) );
				external_function_param_casadi *expl_ode_fun = (external_function_param_casadi *) ptr[0];
				// expl_vde_for
				ptr = (long long *) mxGetData( mxGetField( prhs[3], 0, "expl_vde_for" ) );
				external_function_param_casadi *expl_vde_for = (external_function_param_casadi *) ptr[0];
				// expl_vde_adj
				ptr = (long long *) mxGetData( mxGetField( prhs[3], 0, "expl_vde_adj" ) );
				external_function_param_casadi *expl_vde_adj = (external_function_param_casadi *) ptr[0];

				expl_ode_fun->set_param(expl_ode_fun, p);
				expl_vde_for->set_param(expl_vde_for, p);
				}
			else if(!strcmp(method, "irk"))
				{
				// impl_ode_fun
				ptr = (long long *) mxGetData( mxGetField( prhs[3], 0, "impl_ode_fun" ) );
				external_function_param_casadi *impl_ode_fun = (external_function_param_casadi *) ptr[0];
				// impl_ode_fun_jac_x_xdot
				ptr = (long long *) mxGetData( mxGetField( prhs[3], 0, "impl_ode_fun_jac_x_xdot" ) );
				external_function_param_casadi *impl_ode_fun_jac_x_xdot = (external_function_param_casadi *) ptr[0];
				// impl_ode_jac_x_xdot_u
				ptr = (long long *) mxGetData( mxGetField( prhs[3], 0, "impl_ode_jac_x_xdot_u" ) );
				external_function_param_casadi *impl_ode_jac_x_xdot_u = (external_function_param_casadi *) ptr[0];

				impl_ode_fun->set_param(impl_ode_fun, p);
				impl_ode_fun_jac_x_xdot->set_param(impl_ode_fun_jac_x_xdot, p);
				impl_ode_jac_x_xdot_u->set_param(impl_ode_jac_x_xdot_u, p);
				}
			else
				{
				mexPrintf("\nsim_set: method not supported %s\n", method);
				return;
				}
			}
		else
			{
			mexPrintf("\nsim_set: can not set p for non-param models\n");
			return;
			}
		}
	else
		{
		mexPrintf("\nsim_set: field not supported: %s\n", field);
		return;
		}



	/* return */

	return;

	}


