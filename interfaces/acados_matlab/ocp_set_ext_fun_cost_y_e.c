// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados_c/ocp_nlp_interface.h"
#include "acados/utils/external_function_generic.h"
#include "acados_c/external_function_interface.h"
// mex
#include "mex.h"



// casadi functions for the model
#include "ocp_model.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin sim_impl_ext_fun_create\n");

	// sizeof(long long) == sizeof(void *) = 64 !!!
	long long *ptr;



	/* RHS */

	// C_ocp

	// C_ocp_ext_fun

	// model

	char *cost_type_e;
	if (mxGetField( prhs[2], 0, "cost_type_e" )!=NULL)
		{
		cost_type_e = mxArrayToString( mxGetField( prhs[2], 0, "cost_type_e" ) );
		}
	int np = 0;
	if(mxGetField( prhs[2], 0, "dim_np" )!=NULL) // TODO bool
		{
		np = mxGetScalar( mxGetField( prhs[2], 0, "dim_np" ) );
		}
	// TODO bool instead !!!
	char *param_y_e = mxArrayToString( mxGetField( prhs[2], 0, "cost_param_y_e" ) );


	// opts

//	bool sens_forw = mxGetScalar( mxGetField( prhs[3], 0, "sens_forw" ) );
//	mexPrintf("\n%d\n", sens_forw);
//	char *sim_method = mxArrayToString( mxGetField( prhs[3], 0, "sim_method" ) );
//	mexPrintf("\n%s\n", sim_method);



	/* LHS */

	/* copy existing fields */

//	plhs[0] = mxCreateSharedDataCopy(prhs[1]);
	plhs[0] = mxDuplicateArray(prhs[1]);



	/* populate new fields */

	external_function_casadi *ext_fun_ptr;
	external_function_param_casadi *ext_fun_param_ptr;

	// TODO templetize the casadi function names !!!
	if(!strcmp(cost_type_e, "nonlinear_ls"))
		{
		if(!strcmp(param_y_e, "true")) // TODO bool
			{
			// y_e_fun_jac_ut_xt
			ext_fun_param_ptr = (external_function_param_casadi *) malloc(1*sizeof(external_function_param_casadi));
			external_function_param_casadi_set_fun(ext_fun_param_ptr, &ocp_model_cost_y_e_fun_jac_ut_xt);
			external_function_param_casadi_set_work(ext_fun_param_ptr, &ocp_model_cost_y_e_fun_jac_ut_xt_work);
			external_function_param_casadi_set_sparsity_in(ext_fun_param_ptr, &ocp_model_cost_y_e_fun_jac_ut_xt_sparsity_in);
			external_function_param_casadi_set_sparsity_out(ext_fun_param_ptr, &ocp_model_cost_y_e_fun_jac_ut_xt_sparsity_out);
			external_function_param_casadi_set_n_in(ext_fun_param_ptr, &ocp_model_cost_y_e_fun_jac_ut_xt_n_in);
			external_function_param_casadi_set_n_out(ext_fun_param_ptr, &ocp_model_cost_y_e_fun_jac_ut_xt_n_out);
			external_function_param_casadi_create(ext_fun_param_ptr, np);
			// populate output struct
			mxArray *y_e_fun_jac_ut_xt_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(y_e_fun_jac_ut_xt_mat);
			ptr[0] = (long long) ext_fun_param_ptr;
			mxSetField(plhs[0], 0, "cost_y_e_fun_jac_ut_xt", y_e_fun_jac_ut_xt_mat);
			}
		else
			{
			// y_e_fun_jac_ut_xt
			ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
			external_function_casadi_set_fun(ext_fun_ptr, &ocp_model_cost_y_e_fun_jac_ut_xt);
			external_function_casadi_set_work(ext_fun_ptr, &ocp_model_cost_y_e_fun_jac_ut_xt_work);
			external_function_casadi_set_sparsity_in(ext_fun_ptr, &ocp_model_cost_y_e_fun_jac_ut_xt_sparsity_in);
			external_function_casadi_set_sparsity_out(ext_fun_ptr, &ocp_model_cost_y_e_fun_jac_ut_xt_sparsity_out);
			external_function_casadi_set_n_in(ext_fun_ptr, &ocp_model_cost_y_e_fun_jac_ut_xt_n_in);
			external_function_casadi_set_n_out(ext_fun_ptr, &ocp_model_cost_y_e_fun_jac_ut_xt_n_out);
			external_function_casadi_create(ext_fun_ptr);
			// populate output struct
			mxArray *y_e_fun_jac_ut_xt_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(y_e_fun_jac_ut_xt_mat);
			ptr[0] = (long long) ext_fun_ptr;
			mxSetField(plhs[0], 0, "cost_y_e_fun_jac_ut_xt", y_e_fun_jac_ut_xt_mat);
			}
		}
	
	return;

	}






