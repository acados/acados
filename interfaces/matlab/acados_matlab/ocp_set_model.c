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

	// opts_struct

//	int num_stages = mxGetScalar( mxGetField( prhs[0], 0, "num_stages" ) );
//	int num_steps = mxGetScalar( mxGetField( prhs[0], 0, "num_steps" ) );
//	bool sens_forw = mxGetScalar( mxGetField( prhs[0], 0, "sens_forw" ) );
//	mexPrintf("\n%d\n", sens_forw);
	char *sim_solver = mxArrayToString( mxGetField( prhs[0], 0, "sim_solver" ) );
//	mexPrintf("\n%s\n", sim_solver);


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

	if(!strcmp(sim_solver, "erk"))
		{
//		mexPrintf("\n%s\n", sim_solver);

		// expl_ode_fun
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "expl_ode_fun" ) );
		external_function_casadi *expl_ode_fun = (external_function_casadi *) ptr[0];
		// expl_vde_for
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "expl_vde_for" ) );
		external_function_casadi *expl_vde_for = (external_function_casadi *) ptr[0];
		// expl_vde_adj
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "expl_vde_adj" ) );
		external_function_casadi *expl_vde_adj = (external_function_casadi *) ptr[0];

//		sim_in_set(config, dims, in, "expl_ode_fun", expl_ode_fun);
//		sim_in_set(config, dims, in, "expl_vde_for", expl_vde_for);
		for(ii=0; ii<N; ii++)
			{
//			status = ocp_nlp_dynamics_model_set(config, in, ii, "expl_ode_fun", &expl_ode_fun[ii]); // not needed as sens_forw=1
			status = ocp_nlp_dynamics_model_set(config, in, ii, "expl_vde_for", expl_vde_for);
			}
		}
	else if(!strcmp(sim_solver, "irk"))
		{
//		mexPrintf("\n%s\n", sim_solver);

		// impl_ode_fun
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "impl_ode_fun" ) );
		external_function_casadi *impl_ode_fun = (external_function_casadi *) ptr[0];
		// impl_ode_fun_jac_x_xdot
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "impl_ode_fun_jac_x_xdot" ) );
		external_function_casadi *impl_ode_fun_jac_x_xdot = (external_function_casadi *) ptr[0];
		// impl_ode_jac_x_xdot_u
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "impl_ode_jac_x_xdot_u" ) );
		external_function_casadi *impl_ode_jac_x_xdot_u = (external_function_casadi *) ptr[0];

//		sim_in_set(config, dims, in, "impl_ode_fun", impl_ode_fun);
//		sim_in_set(config, dims, in, "impl_ode_fun_jac_x_xdot", impl_ode_fun_jac_x_xdot);
//		sim_in_set(config, dims, in, "impl_ode_jac_x_xdot_u", impl_ode_jac_x_xdot_u);
		for(ii=0; ii<N; ii++)
			{
			status = ocp_nlp_dynamics_model_set(config, in, ii, "impl_ode_fun", impl_ode_fun);
			status = ocp_nlp_dynamics_model_set(config, in, ii, "impl_ode_fun_jac_x_xdot", impl_ode_fun_jac_x_xdot);
			status = ocp_nlp_dynamics_model_set(config, in, ii, "impl_ode_jac_x_xdot_u", impl_ode_jac_x_xdot_u);
			}
		}
	else
		{
		mexPrintf("\nsim_set_model: sim_solver not supported %s\n", sim_solver);
		return;
		}



	/* return */

	return;

	}


