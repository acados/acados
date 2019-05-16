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



// casadi functions for the model
#include "sim_model.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin sim_expl_ext_fun_create\n");

	// sizeof(long long) == sizeof(void *) = 64 !!!
	long long *ptr;



	/* RHS */

	// C_sim

	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "config" ) );
	sim_config *config = (sim_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "dims" ) );
	void *dims = (void *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "in" ) );
	sim_in *in = (sim_in *) ptr[0];

	// model

	int np = 0;
	if(mxGetField( prhs[2], 0, "dim_np" )!=NULL) // TODO bool
		{
		np = mxGetScalar( mxGetField( prhs[2], 0, "dim_np" ) );
		}
	// TODO bool instead !!!
	char *param_f = mxArrayToString( mxGetField( prhs[2], 0, "dyn_param_f" ) );

	// opts

//	bool sens_forw = mxGetScalar( mxGetField( prhs[3], 0, "sens_forw" ) );
//	mexPrintf("\n%d\n", sens_forw);
	char *method = mxArrayToString( mxGetField( prhs[3], 0, "method" ) );
//	mexPrintf("\n%s\n", method);
//	bool sens_adj = mxGetScalar( mxGetField( prhs[3], 0, "sens_adj" ) );
//	mexPrintf("\n%d\n", sens_adj);



	/* LHS */

	/* copy existing fields */

//	plhs[0] = mxCreateSharedDataCopy(prhs[1]);
	plhs[0] = mxDuplicateArray(prhs[1]);



	/* populate new fields */

	external_function_casadi *ext_fun_ptr;
	external_function_param_casadi *ext_fun_param_ptr;

	mxArray *tmp_mat;

	// TODO templetize the casadi function names !!!

	// TODO set checking opts !!!

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
			sim_in_set(config, dims, in, "expl_ode_fun", ext_fun_param_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_param_ptr;
			mxSetField(plhs[0], 0, "dyn_expl_ode_fun", tmp_mat);

			// expl_vde_for
			ext_fun_param_ptr = (external_function_param_casadi *) malloc(1*sizeof(external_function_param_casadi));
			external_function_param_casadi_set_fun(ext_fun_param_ptr, &sim_model_dyn_expl_vde_for);
			external_function_param_casadi_set_work(ext_fun_param_ptr, &sim_model_dyn_expl_vde_for_work);
			external_function_param_casadi_set_sparsity_in(ext_fun_param_ptr, &sim_model_dyn_expl_vde_for_sparsity_in);
			external_function_param_casadi_set_sparsity_out(ext_fun_param_ptr, &sim_model_dyn_expl_vde_for_sparsity_out);
			external_function_param_casadi_set_n_in(ext_fun_param_ptr, &sim_model_dyn_expl_vde_for_n_in);
			external_function_param_casadi_set_n_out(ext_fun_param_ptr, &sim_model_dyn_expl_vde_for_n_out);
			external_function_param_casadi_create(ext_fun_param_ptr, np);
			sim_in_set(config, dims, in, "expl_vde_for", ext_fun_param_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_param_ptr;
			mxSetField(plhs[0], 0, "dyn_expl_vde_for", tmp_mat);

			// expl_vde_adj
			ext_fun_param_ptr = (external_function_param_casadi *) malloc(1*sizeof(external_function_param_casadi));
			external_function_param_casadi_set_fun(ext_fun_param_ptr, &sim_model_dyn_expl_vde_adj);
			external_function_param_casadi_set_work(ext_fun_param_ptr, &sim_model_dyn_expl_vde_adj_work);
			external_function_param_casadi_set_sparsity_in(ext_fun_param_ptr, &sim_model_dyn_expl_vde_adj_sparsity_in);
			external_function_param_casadi_set_sparsity_out(ext_fun_param_ptr, &sim_model_dyn_expl_vde_adj_sparsity_out);
			external_function_param_casadi_set_n_in(ext_fun_param_ptr, &sim_model_dyn_expl_vde_adj_n_in);
			external_function_param_casadi_set_n_out(ext_fun_param_ptr, &sim_model_dyn_expl_vde_adj_n_out);
			external_function_param_casadi_create(ext_fun_param_ptr, np);
			sim_in_set(config, dims, in, "expl_vde_adj", ext_fun_param_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_param_ptr;
			mxSetField(plhs[0], 0, "dyn_expl_vde_adj", tmp_mat);

			// expl_ode_hes
			ext_fun_param_ptr = (external_function_param_casadi *) malloc(1*sizeof(external_function_param_casadi));
			external_function_param_casadi_set_fun(ext_fun_param_ptr, &sim_model_dyn_expl_ode_hes);
			external_function_param_casadi_set_work(ext_fun_param_ptr, &sim_model_dyn_expl_ode_hes_work);
			external_function_param_casadi_set_sparsity_in(ext_fun_param_ptr, &sim_model_dyn_expl_ode_hes_sparsity_in);
			external_function_param_casadi_set_sparsity_out(ext_fun_param_ptr, &sim_model_dyn_expl_ode_hes_sparsity_out);
			external_function_param_casadi_set_n_in(ext_fun_param_ptr, &sim_model_dyn_expl_ode_hes_n_in);
			external_function_param_casadi_set_n_out(ext_fun_param_ptr, &sim_model_dyn_expl_ode_hes_n_out);
			external_function_param_casadi_create(ext_fun_param_ptr, np);
			sim_in_set(config, dims, in, "expl_ode_hes", ext_fun_param_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_param_ptr;
			mxSetField(plhs[0], 0, "dyn_expl_ode_hes", tmp_mat);
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
			sim_in_set(config, dims, in, "expl_ode_fun", ext_fun_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_ptr;
			mxSetField(plhs[0], 0, "dyn_expl_ode_fun", tmp_mat);

			// expl_vde_for
			ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
			external_function_casadi_set_fun(ext_fun_ptr, &sim_model_dyn_expl_vde_for);
			external_function_casadi_set_work(ext_fun_ptr, &sim_model_dyn_expl_vde_for_work);
			external_function_casadi_set_sparsity_in(ext_fun_ptr, &sim_model_dyn_expl_vde_for_sparsity_in);
			external_function_casadi_set_sparsity_out(ext_fun_ptr, &sim_model_dyn_expl_vde_for_sparsity_out);
			external_function_casadi_set_n_in(ext_fun_ptr, &sim_model_dyn_expl_vde_for_n_in);
			external_function_casadi_set_n_out(ext_fun_ptr, &sim_model_dyn_expl_vde_for_n_out);
			external_function_casadi_create(ext_fun_ptr);
			sim_in_set(config, dims, in, "expl_vde_for", ext_fun_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_ptr;
			mxSetField(plhs[0], 0, "dyn_expl_vde_for", tmp_mat);

			// expl_vde_adj
			ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
			external_function_casadi_set_fun(ext_fun_ptr, &sim_model_dyn_expl_vde_adj);
			external_function_casadi_set_work(ext_fun_ptr, &sim_model_dyn_expl_vde_adj_work);
			external_function_casadi_set_sparsity_in(ext_fun_ptr, &sim_model_dyn_expl_vde_adj_sparsity_in);
			external_function_casadi_set_sparsity_out(ext_fun_ptr, &sim_model_dyn_expl_vde_adj_sparsity_out);
			external_function_casadi_set_n_in(ext_fun_ptr, &sim_model_dyn_expl_vde_adj_n_in);
			external_function_casadi_set_n_out(ext_fun_ptr, &sim_model_dyn_expl_vde_adj_n_out);
			external_function_casadi_create(ext_fun_ptr);
			sim_in_set(config, dims, in, "expl_vde_adj", ext_fun_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_ptr;
			mxSetField(plhs[0], 0, "dyn_expl_vde_adj", tmp_mat);

			// expl_ode_hes
			ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
			external_function_casadi_set_fun(ext_fun_ptr, &sim_model_dyn_expl_ode_hes);
			external_function_casadi_set_work(ext_fun_ptr, &sim_model_dyn_expl_ode_hes_work);
			external_function_casadi_set_sparsity_in(ext_fun_ptr, &sim_model_dyn_expl_ode_hes_sparsity_in);
			external_function_casadi_set_sparsity_out(ext_fun_ptr, &sim_model_dyn_expl_ode_hes_sparsity_out);
			external_function_casadi_set_n_in(ext_fun_ptr, &sim_model_dyn_expl_ode_hes_n_in);
			external_function_casadi_set_n_out(ext_fun_ptr, &sim_model_dyn_expl_ode_hes_n_out);
			external_function_casadi_create(ext_fun_ptr);
			sim_in_set(config, dims, in, "expl_ode_hes", ext_fun_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_ptr;
			mxSetField(plhs[0], 0, "dyn_expl_ode_hes", tmp_mat);
			}
		}
	else
		{
		mexPrintf("\nsim_set_ext_fun_dyn_expl: method not supported %s\n", method);
		return;
		}
	
	return;

	}



