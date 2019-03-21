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


	// model

	// TODO bool instead !!!
	char *param_f = mxArrayToString( mxGetField( prhs[0], 0, "dyn_param_f" ) );
	char *param_h = mxArrayToString( mxGetField( prhs[0], 0, "constr_param_h" ) );
	char *param_h_e = mxArrayToString( mxGetField( prhs[0], 0, "constr_param_h_e" ) );
	char *param_y = mxArrayToString( mxGetField( prhs[0], 0, "cost_param_y" ) );
	char *param_y_e = mxArrayToString( mxGetField( prhs[0], 0, "cost_param_y_e" ) );
	char *param_ext_cost = mxArrayToString( mxGetField( prhs[0], 0, "cost_param_ext_cost" ) );
	char *param_ext_cost_e = mxArrayToString( mxGetField( prhs[0], 0, "cost_param_ext_cost_e" ) );

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



	int N = dims->N;


	// C_sim_ext_fun

	
	// TODO check for empty struct member


	/* set in model */

	external_function_casadi *ext_fun_ptr;
	external_function_param_casadi *ext_fun_param_ptr;

	if (mxGetField( prhs[1], 0, "dyn_expl_ode_fun" )!=NULL && mxGetM(mxGetField( prhs[1], 0, "dyn_expl_ode_fun" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dyn_expl_ode_fun" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "expl_ode_fun", ext_fun_param_ptr+ii); // NOTE not needed as sens_forw=1
				}
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "expl_ode_fun", ext_fun_ptr+ii); // NOTE not needed as sens_forw=1
				}
			}
		}
	if (mxGetField( prhs[1], 0, "dyn_expl_vde_for" )!=NULL && mxGetM(mxGetField( prhs[1], 0, "dyn_expl_vde_for" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dyn_expl_vde_for" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "expl_vde_for", ext_fun_param_ptr+ii);
				}
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "expl_vde_for", ext_fun_ptr+ii);
				}
			}
		}
	if (mxGetField( prhs[1], 0, "dyn_expl_vde_adj" )!=NULL && mxGetM(mxGetField( prhs[1], 0, "dyn_expl_vde_adj" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dyn_expl_vde_adj" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "expl_vde_adj", ext_fun_param_ptr+ii);
				}
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "expl_vde_adj", ext_fun_ptr+ii);
				}
			}
		}
	if (mxGetField( prhs[1], 0, "dyn_impl_ode_fun" )!=NULL && mxGetM(mxGetField( prhs[1], 0, "dyn_impl_ode_fun" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dyn_impl_ode_fun" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "impl_ode_fun", ext_fun_param_ptr+ii);
				}
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "impl_ode_fun", ext_fun_ptr+ii);
				}
			}
		}
	if (mxGetField( prhs[1], 0, "dyn_impl_ode_fun_jac_x_xdot" )!=NULL && mxGetM(mxGetField( prhs[1], 0, "dyn_impl_ode_fun_jac_x_xdot" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dyn_impl_ode_fun_jac_x_xdot" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "impl_ode_fun_jac_x_xdot", ext_fun_param_ptr+ii);
				}
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "impl_ode_fun_jac_x_xdot", ext_fun_ptr+ii);
				}
			}
		}
	if (mxGetField( prhs[1], 0, "dyn_impl_ode_jac_x_xdot_u" )!=NULL && mxGetM(mxGetField( prhs[1], 0, "dyn_impl_ode_jac_x_xdot_u" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dyn_impl_ode_jac_x_xdot_u" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "impl_ode_jac_x_xdot_u", ext_fun_param_ptr+ii);
				}
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_dynamics_model_set(config, dims, in, ii, "impl_ode_jac_x_xdot_u", ext_fun_ptr+ii);
				}
			}
		}
	if (mxGetField( prhs[1], 0, "constr_h_fun_jac_ut_xt" )!=NULL && mxGetM(mxGetField( prhs[1], 0, "constr_h_fun_jac_ut_xt" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "constr_h_fun_jac_ut_xt" ) );
		if(!strcmp(param_h, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_constraints_model_set(config, dims, in, ii, "nl_constr_h_fun_jac", ext_fun_param_ptr+ii);
				}
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_constraints_model_set(config, dims, in, ii, "nl_constr_h_fun_jac", ext_fun_ptr+ii);
				}
			}
		}
	if (mxGetField( prhs[1], 0, "constr_h_e_fun_jac_ut_xt" )!=NULL && mxGetM(mxGetField( prhs[1], 0, "constr_h_e_fun_jac_ut_xt" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "constr_h_e_fun_jac_ut_xt" ) );
		if(!strcmp(param_h_e, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			status = ocp_nlp_constraints_model_set(config, dims, in, N, "nl_constr_h_fun_jac", ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			status = ocp_nlp_constraints_model_set(config, dims, in, N, "nl_constr_h_fun_jac", ext_fun_ptr);
			}
		}
	if (mxGetField( prhs[1], 0, "cost_y_fun_jac_ut_xt" )!=NULL && mxGetM(mxGetField( prhs[1], 0, "cost_y_fun_jac_ut_xt" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "cost_y_fun_jac_ut_xt" ) );
		if(!strcmp(param_y, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_cost_model_set(config, dims, in, ii, "nls_res_jac", ext_fun_param_ptr+ii);
				}
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_cost_model_set(config, dims, in, ii, "nls_res_jac", ext_fun_ptr+ii);
				}
			}
		}
	if (mxGetField( prhs[1], 0, "cost_y_e_fun_jac_ut_xt" )!=NULL && mxGetM(mxGetField( prhs[1], 0, "cost_y_e_fun_jac_ut_xt" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "cost_y_e_fun_jac_ut_xt" ) );
		if(!strcmp(param_y_e, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			status = ocp_nlp_cost_model_set(config, dims, in, N, "nls_res_jac", ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			status = ocp_nlp_cost_model_set(config, dims, in, N, "nls_res_jac", ext_fun_ptr);
			}
		}
	if (mxGetField( prhs[1], 0, "cost_ext_cost_jac_hes" )!=NULL && mxGetM(mxGetField( prhs[1], 0, "cost_ext_cost_jac_hes" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "cost_ext_cost_jac_hes" ) );
		if(!strcmp(param_ext_cost, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_cost_model_set(config, dims, in, ii, "ext_cost_jac_hes", ext_fun_param_ptr+ii);
				}
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				status = ocp_nlp_cost_model_set(config, dims, in, ii, "ext_cost_jac_hes", ext_fun_ptr+ii);
				}
			}
		}
	if (mxGetField( prhs[1], 0, "cost_ext_cost_e_jac_hes" )!=NULL && mxGetM(mxGetField( prhs[1], 0, "cost_ext_cost_e_jac_hes" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "cost_ext_cost_e_jac_hes" ) );
		if(!strcmp(param_ext_cost_e, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			status = ocp_nlp_cost_model_set(config, dims, in, N, "ext_cost_jac_hes", ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			status = ocp_nlp_cost_model_set(config, dims, in, N, "ext_cost_jac_hes", ext_fun_ptr);
			}
		}



	/* return */

	return;

	}


