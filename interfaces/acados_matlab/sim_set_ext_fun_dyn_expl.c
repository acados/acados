// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados/utils/external_function_generic.h"
#include "acados_c/external_function_interface.h"
// mex
#include "mex.h"



// casadi functions for the model
#include "sim_model.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin sim_expl_ext_fun_create\n");

	// sizeof(long long) == sizeof(void *) = 64 !!!
	long long *ptr;



	/* RHS */

	// model

	int np = 0;
	if(mxGetField( prhs[1], 0, "dim_np" )!=NULL) // TODO bool
		{
		np = mxGetScalar( mxGetField( prhs[1], 0, "dim_np" ) );
		}
	// TODO bool instead !!!
	char *param_f = mxArrayToString( mxGetField( prhs[1], 0, "dyn_param_f" ) );

	// opts

	// TODO use them !!!
//	bool sens_forw = mxGetScalar( mxGetField( prhs[0], 0, "sens_forw" ) );
//	mexPrintf("\n%d\n", sens_forw);
	char *method = mxArrayToString( mxGetField( prhs[2], 0, "method" ) );
//	mexPrintf("\n%s\n", method);



	/* LHS */

	/* copy existing fields */

//	plhs[0] = mxCreateSharedDataCopy(prhs[0]);
	plhs[0] = mxDuplicateArray(prhs[0]);



	/* populate new fields */

	external_function_casadi *ext_fun_ptr;
	external_function_param_casadi *ext_fun_param_ptr;

	// TODO templetize the casadi function names !!!
	if(!strcmp(method, "erk"))
		{
		if(!strcmp(param_f, "true")) // TODO bool
			{
			// expl_ode_fun
			ext_fun_param_ptr = (external_function_param_casadi *) malloc(1*sizeof(external_function_param_casadi));
			external_function_param_casadi_set_fun(ext_fun_param_ptr, &sim_model_dyn_expl_ode_fun);
			external_function_param_casadi_set_work(ext_fun_param_ptr, &sim_model_dyn_expl_ode_fun_work);
			external_function_param_casadi_set_sparsity_in(ext_fun_param_ptr, &sim_model_dyn_expl_ode_fun_sparsity_in);
			external_function_param_casadi_set_sparsity_out(ext_fun_param_ptr, &sim_model_dyn_expl_ode_fun_sparsity_out);
			external_function_param_casadi_set_n_in(ext_fun_param_ptr, &sim_model_dyn_expl_ode_fun_n_in);
			external_function_param_casadi_set_n_out(ext_fun_param_ptr, &sim_model_dyn_expl_ode_fun_n_out);
			external_function_param_casadi_create(ext_fun_param_ptr, np);
			// populate output struct
			mxArray *expl_ode_fun_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(expl_ode_fun_mat);
			ptr[0] = (long long) ext_fun_param_ptr;
			mxSetField(plhs[0], 0, "dyn_expl_ode_fun", expl_ode_fun_mat);

			// expl_vde_for
			ext_fun_param_ptr = (external_function_param_casadi *) malloc(1*sizeof(external_function_param_casadi));
			external_function_param_casadi_set_fun(ext_fun_param_ptr, &sim_model_dyn_expl_vde_for);
			external_function_param_casadi_set_work(ext_fun_param_ptr, &sim_model_dyn_expl_vde_for_work);
			external_function_param_casadi_set_sparsity_in(ext_fun_param_ptr, &sim_model_dyn_expl_vde_for_sparsity_in);
			external_function_param_casadi_set_sparsity_out(ext_fun_param_ptr, &sim_model_dyn_expl_vde_for_sparsity_out);
			external_function_param_casadi_set_n_in(ext_fun_param_ptr, &sim_model_dyn_expl_vde_for_n_in);
			external_function_param_casadi_set_n_out(ext_fun_param_ptr, &sim_model_dyn_expl_vde_for_n_out);
			external_function_param_casadi_create(ext_fun_param_ptr, np);
			// populate output struct
			mxArray *expl_vde_for_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(expl_vde_for_mat);
			ptr[0] = (long long) ext_fun_param_ptr;
			mxSetField(plhs[0], 0, "dyn_expl_vde_for", expl_vde_for_mat);
			}
		else
			{
			// expl_ode_fun
			ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
			external_function_casadi_set_fun(ext_fun_ptr, &sim_model_dyn_expl_ode_fun);
			external_function_casadi_set_work(ext_fun_ptr, &sim_model_dyn_expl_ode_fun_work);
			external_function_casadi_set_sparsity_in(ext_fun_ptr, &sim_model_dyn_expl_ode_fun_sparsity_in);
			external_function_casadi_set_sparsity_out(ext_fun_ptr, &sim_model_dyn_expl_ode_fun_sparsity_out);
			external_function_casadi_set_n_in(ext_fun_ptr, &sim_model_dyn_expl_ode_fun_n_in);
			external_function_casadi_set_n_out(ext_fun_ptr, &sim_model_dyn_expl_ode_fun_n_out);
			external_function_casadi_create(ext_fun_ptr);
			// populate output struct
			mxArray *expl_ode_fun_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(expl_ode_fun_mat);
			ptr[0] = (long long) ext_fun_ptr;
			mxSetField(plhs[0], 0, "dyn_expl_ode_fun", expl_ode_fun_mat);

			// expl_vde_for
			ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
			external_function_casadi_set_fun(ext_fun_ptr, &sim_model_dyn_expl_vde_for);
			external_function_casadi_set_work(ext_fun_ptr, &sim_model_dyn_expl_vde_for_work);
			external_function_casadi_set_sparsity_in(ext_fun_ptr, &sim_model_dyn_expl_vde_for_sparsity_in);
			external_function_casadi_set_sparsity_out(ext_fun_ptr, &sim_model_dyn_expl_vde_for_sparsity_out);
			external_function_casadi_set_n_in(ext_fun_ptr, &sim_model_dyn_expl_vde_for_n_in);
			external_function_casadi_set_n_out(ext_fun_ptr, &sim_model_dyn_expl_vde_for_n_out);
			external_function_casadi_create(ext_fun_ptr);
			// populate output struct
			mxArray *expl_vde_for_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(expl_vde_for_mat);
			ptr[0] = (long long) ext_fun_ptr;
			mxSetField(plhs[0], 0, "dyn_expl_vde_for", expl_vde_for_mat);
			}
		}
	else
		{
		mexPrintf("\nsim_set_ext_fun_dyn_expl: method not supported %s\n", method);
		return;
		}
mexPrintf("\nhere\n");
	
	return;

	}



