// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
//#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados/utils/external_function_generic.h"
#include "acados_c/external_function_interface.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin sim_set_model\n");

	long long *ptr;

	int ii, status;


	/* RHS */


	// C_ocp

	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "config" ) );
	ocp_nlp_config *config = (ocp_nlp_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dims" ) );
	ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "in" ) );
	ocp_nlp_in *in = (ocp_nlp_in *) ptr[0];



	int N = dims->N;


	// C_sim_ext_fun

	
	// TODO check for empty struct member


	/* set in model */

	external_function_casadi *ext_fun_ptr;

	if (mxGetField( prhs[0], 0, "expl_ode_fun" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "expl_ode_fun" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		for(ii=0; ii<N; ii++)
			{
			status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "expl_ode_fun", ext_fun_ptr); // NOTE not needed as sens_forw=1
			}
		}
	if (mxGetField( prhs[0], 0, "expl_vde_for" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "expl_vde_for" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		for(ii=0; ii<N; ii++)
			{
			status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "expl_vde_for", ext_fun_ptr);
			}
		}
	if (mxGetField( prhs[0], 0, "expl_vde_adj" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "expl_vde_adj" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		for(ii=0; ii<N; ii++)
			{
			status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "expl_vde_adj", ext_fun_ptr);
			}
		}
	if (mxGetField( prhs[0], 0, "impl_ode_fun" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "impl_ode_fun" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		for(ii=0; ii<N; ii++)
			{
			status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "impl_ode_fun", ext_fun_ptr);
			}
		}
	if (mxGetField( prhs[0], 0, "impl_ode_fun_jac_x_xdot" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "impl_ode_fun_jac_x_xdot" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		for(ii=0; ii<N; ii++)
			{
			status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "impl_ode_fun_jac_x_xdot", ext_fun_ptr);
			}
		}
	if (mxGetField( prhs[0], 0, "impl_ode_jac_x_xdot_u" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "impl_ode_jac_x_xdot_u" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		for(ii=0; ii<N; ii++)
			{
			status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "impl_ode_jac_x_xdot_u", ext_fun_ptr);
			}
		}
	if (mxGetField( prhs[0], 0, "h_fun_jac_ut_xt" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "h_fun_jac_ut_xt" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		for(ii=0; ii<N; ii++)
			{
			status = ocp_nlp_constraints_model_set(config, dims, in, ii, "h", ext_fun_ptr);
			}
		}
	if (mxGetField( prhs[0], 0, "h_e_fun_jac_ut_xt" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "h_e_fun_jac_ut_xt" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		status = ocp_nlp_constraints_model_set(config, dims, in, N, "h", ext_fun_ptr);
		}
	if (mxGetField( prhs[0], 0, "y_fun_jac_ut_xt" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "y_fun_jac_ut_xt" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		for(ii=0; ii<N; ii++)
			{
			status = ocp_nlp_cost_model_set(config, dims, in, ii, "nls_res_jac", ext_fun_ptr);
			}
		}
	if (mxGetField( prhs[0], 0, "y_e_fun_jac_ut_xt" )!=NULL)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "y_e_fun_jac_ut_xt" ) );
		ext_fun_ptr = (external_function_casadi *) ptr[0];
		status = ocp_nlp_cost_model_set(config, dims, in, N, "nls_res_jac", ext_fun_ptr);
		}



	/* return */

	return;

	}


