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
	int status;
	long long *ptr;



	/* RHS */

	// C_ocp

	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "config" ) );
	ocp_nlp_config *config = (ocp_nlp_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "dims" ) );
	ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "in" ) );
	ocp_nlp_in *in = (ocp_nlp_in *) ptr[0];

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
	char *param_ext_cost_e = mxArrayToString( mxGetField( prhs[2], 0, "cost_param_ext_cost_e" ) );


	// opts

//	bool sens_forw = mxGetScalar( mxGetField( prhs[3], 0, "sens_forw" ) );
//	mexPrintf("\n%d\n", sens_forw);
//	char *sim_method = mxArrayToString( mxGetField( prhs[3], 0, "sim_method" ) );
//	mexPrintf("\n%s\n", sim_method);


	int N = dims->N;



	/* LHS */

	/* copy existing fields */

//	plhs[0] = mxCreateSharedDataCopy(prhs[1]);
	plhs[0] = mxDuplicateArray(prhs[1]);



	/* populate new fields */

	external_function_casadi *ext_fun_ptr;
	external_function_param_casadi *ext_fun_param_ptr;

	mxArray *tmp_mat;

	// TODO templetize the casadi function names !!!
	if(!strcmp(cost_type_e, "ext_cost"))
		{
		if(!strcmp(param_ext_cost_e, "true")) // TODO bool
			{
			// ext_cost_e_jac_hes
			ext_fun_param_ptr = (external_function_param_casadi *) malloc(1*sizeof(external_function_param_casadi));
			external_function_param_casadi_set_fun(ext_fun_param_ptr, &ocp_model_cost_ext_cost_e_jac_hes);
			external_function_param_casadi_set_work(ext_fun_param_ptr, &ocp_model_cost_ext_cost_e_jac_hes_work);
			external_function_param_casadi_set_sparsity_in(ext_fun_param_ptr, &ocp_model_cost_ext_cost_e_jac_hes_sparsity_in);
			external_function_param_casadi_set_sparsity_out(ext_fun_param_ptr, &ocp_model_cost_ext_cost_e_jac_hes_sparsity_out);
			external_function_param_casadi_set_n_in(ext_fun_param_ptr, &ocp_model_cost_ext_cost_e_jac_hes_n_in);
			external_function_param_casadi_set_n_out(ext_fun_param_ptr, &ocp_model_cost_ext_cost_e_jac_hes_n_out);
			external_function_param_casadi_create(ext_fun_param_ptr, np);
			status = ocp_nlp_cost_model_set(config, dims, in, N, "ext_cost_jac_hes", ext_fun_param_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_param_ptr;
			mxSetField(plhs[0], 0, "cost_ext_cost_e_jac_hes", tmp_mat);
			}
		else
			{
			// ext_cost_e_jac_hes
			ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
			external_function_casadi_set_fun(ext_fun_ptr, &ocp_model_cost_ext_cost_e_jac_hes);
			external_function_casadi_set_work(ext_fun_ptr, &ocp_model_cost_ext_cost_e_jac_hes_work);
			external_function_casadi_set_sparsity_in(ext_fun_ptr, &ocp_model_cost_ext_cost_e_jac_hes_sparsity_in);
			external_function_casadi_set_sparsity_out(ext_fun_ptr, &ocp_model_cost_ext_cost_e_jac_hes_sparsity_out);
			external_function_casadi_set_n_in(ext_fun_ptr, &ocp_model_cost_ext_cost_e_jac_hes_n_in);
			external_function_casadi_set_n_out(ext_fun_ptr, &ocp_model_cost_ext_cost_e_jac_hes_n_out);
			external_function_casadi_create(ext_fun_ptr);
			status = ocp_nlp_cost_model_set(config, dims, in, N, "ext_cost_jac_hes", ext_fun_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_ptr;
			mxSetField(plhs[0], 0, "cost_ext_cost_e_jac_hes", tmp_mat);
			}
		}
	
	return;

	}







