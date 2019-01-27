// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados_c/ocp_nlp_interface.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin ocp_get\n");

	long long *ptr;

	int ii;

	/* RHS */

	// model
	// TODO bool instead !!!
	char *param_f = mxArrayToString( mxGetField( prhs[0], 0, "param_f" ) );
	char *param_y = mxArrayToString( mxGetField( prhs[0], 0, "param_y" ) );
	char *param_y_e = mxArrayToString( mxGetField( prhs[0], 0, "param_y_e" ) );
	char *param_h = mxArrayToString( mxGetField( prhs[0], 0, "param_h" ) );
	char *param_h_e = mxArrayToString( mxGetField( prhs[0], 0, "param_h_e" ) );

	// opts
	// TODO
	char *sim_method = mxArrayToString( mxGetField( prhs[1], 0, "sim_method" ) );

	// C_ocp

	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "config" ) );
	ocp_nlp_config *config = (ocp_nlp_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dims" ) );
	ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "in" ) );
	ocp_nlp_in *in = (ocp_nlp_in *) ptr[0];
	// out
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "out" ) );
	ocp_nlp_out *out = (ocp_nlp_out *) ptr[0];

	// field
	char *field = mxArrayToString( prhs[4] );
//	mexPrintf("\n%s\n", field);



	int N = dims->N;
	int nu = dims->nu[0];
	int nx = dims->nx[0];



	// TODO implement with LHS !!!!!
	// value
	if (!strcmp(field, "x0"))
		{
		double *x0 = mxGetPr( prhs[5] );
		ocp_nlp_constraints_model_set(config, dims, in, 0, "lbx", x0);
		ocp_nlp_constraints_model_set(config, dims, in, 0, "ubx", x0);
		}
	else if (!strcmp(field, "yr"))
		{
		double *yr = mxGetPr( prhs[5] );
		for (ii=0; ii<N; ii++)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "y_ref", yr);
			}
		}
	else if (!strcmp(field, "yr_e"))
		{
		double *yr_e = mxGetPr( prhs[5] );
		ocp_nlp_cost_model_set(config, dims, in, N, "y_ref", yr_e);
		}
	else if (!strcmp(field, "x_init"))
		{
		double *x_init = mxGetPr( prhs[5] );
		for (ii=0; ii<=N; ii++)
			{
			ocp_nlp_out_set(config, dims, out, ii, "x", x_init+ii*nx);
			}
		}
	else if (!strcmp(field, "u_init"))
		{
		double *u_init = mxGetPr( prhs[5] );
		for (ii=0; ii<N; ii++)
			{
			ocp_nlp_out_set(config, dims, out, ii, "u", u_init+ii*nu);
			}
		}
	// TODO make setters for ALL numerical data
	else if(!strcmp(field, "p"))
		{
		double *p = mxGetPr( prhs[5] );
		external_function_param_casadi *ext_fun_param_ptr;
		if(!strcmp(param_f, "true")) // TODO bool
			{
			if(!strcmp(sim_method, "erk"))
				{
				// expl_ode_fun
				ptr = (long long *) mxGetData( mxGetField( prhs[3], 0, "expl_ode_fun" ) );
				ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
				ext_fun_param_ptr->set_param(ext_fun_param_ptr, p);
				// expl_vde_for
				ptr = (long long *) mxGetData( mxGetField( prhs[3], 0, "expl_vde_for" ) );
				ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
				ext_fun_param_ptr->set_param(ext_fun_param_ptr, p);
				}
			else if(!strcmp(sim_method, "irk"))
				{
				// impl_ode_fun
				ptr = (long long *) mxGetData( mxGetField( prhs[3], 0, "impl_ode_fun" ) );
				ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
				ext_fun_param_ptr->set_param(ext_fun_param_ptr, p);
				// impl_ode_fun_jac_x_xdot
				ptr = (long long *) mxGetData( mxGetField( prhs[3], 0, "impl_ode_fun_jac_x_xdot" ) );
				ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
				ext_fun_param_ptr->set_param(ext_fun_param_ptr, p);
				// impl_ode_jac_x_xdot_u
				ptr = (long long *) mxGetData( mxGetField( prhs[3], 0, "impl_ode_jac_x_xdot_u" ) );
				ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
				ext_fun_param_ptr->set_param(ext_fun_param_ptr, p);
				}
			else
				{
				mexPrintf("\nocp_set: sim_method not supported %s\n", sim_method);
				return;
				}
			}
		if(!strcmp(param_h, "true")) // TODO bool
			{
			// h_fun_jac_ut_xt
			ptr = (long long *) mxGetData( mxGetField( prhs[3], 0, "h_fun_jac_ut_xt" ) );
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			ext_fun_param_ptr->set_param(ext_fun_param_ptr, p);
			}
		if(!strcmp(param_h_e, "true")) // TODO bool
			{
			// h_e_fun_jac_ut_xt
			ptr = (long long *) mxGetData( mxGetField( prhs[3], 0, "h_e_fun_jac_ut_xt" ) );
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			ext_fun_param_ptr->set_param(ext_fun_param_ptr, p);
			}
		if(!strcmp(param_y, "true")) // TODO bool
			{
			// y_fun_jac_ut_xt
			ptr = (long long *) mxGetData( mxGetField( prhs[3], 0, "y_fun_jac_ut_xt" ) );
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			ext_fun_param_ptr->set_param(ext_fun_param_ptr, p);
			}
		if(!strcmp(param_y_e, "true")) // TODO bool
			{
			// y_e_fun_jac_ut_xt
			ptr = (long long *) mxGetData( mxGetField( prhs[3], 0, "y_e_fun_jac_ut_xt" ) );
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			ext_fun_param_ptr->set_param(ext_fun_param_ptr, p);
			}
//		else
//			{
//			mexPrintf("\nocp_set: can not set p for non-param f, y, y_e, h or h_e\n");
//			return;
//			}
		}
	else
		{
		mexPrintf("\nocp_set: field not supported: %s\n", field);
		return;
		}



	/* return */

	return;

	}




