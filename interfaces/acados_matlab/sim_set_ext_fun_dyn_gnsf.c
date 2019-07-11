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

//	mexPrintf("\nin sim_impl_ext_fun_create\n");

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

	if(!strcmp(method, "irk_gnsf"))
		{
		if(!strcmp(param_f, "true")) // TODO bool
			{
			// gnsf_f_lo_fun_jac_x1k1uz
			ext_fun_param_ptr = (external_function_param_casadi *) malloc(1*sizeof(external_function_param_casadi));
			external_function_param_casadi_set_fun(ext_fun_param_ptr, &sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz);
			external_function_param_casadi_set_work(ext_fun_param_ptr, &sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz_work);
			external_function_param_casadi_set_sparsity_in(ext_fun_param_ptr, &sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz_sparsity_in);
			external_function_param_casadi_set_sparsity_out(ext_fun_param_ptr, &sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz_sparsity_out);
			external_function_param_casadi_set_n_in(ext_fun_param_ptr, &sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz_n_in);
			external_function_param_casadi_set_n_out(ext_fun_param_ptr, &sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz_n_out);
			external_function_param_casadi_create(ext_fun_param_ptr, np);
			sim_in_set(config, dims, in, "gnsf_f_lo_fun_jac_x1k1uz", ext_fun_param_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_param_ptr;
			mxSetField(plhs[0], 0, "dyn_gnsf_f_lo_fun_jac_x1k1uz", tmp_mat);

			// gnsf_get_matrices_fun
			ext_fun_param_ptr = (external_function_param_casadi *) malloc(1*sizeof(external_function_param_casadi));
			external_function_param_casadi_set_fun(ext_fun_param_ptr, &sim_model_dyn_gnsf_get_matrices_fun);
			external_function_param_casadi_set_work(ext_fun_param_ptr, &sim_model_dyn_gnsf_get_matrices_fun_work);
			external_function_param_casadi_set_sparsity_in(ext_fun_param_ptr, &sim_model_dyn_gnsf_get_matrices_fun_sparsity_in);
			external_function_param_casadi_set_sparsity_out(ext_fun_param_ptr, &sim_model_dyn_gnsf_get_matrices_fun_sparsity_out);
			external_function_param_casadi_set_n_in(ext_fun_param_ptr, &sim_model_dyn_gnsf_get_matrices_fun_n_in);
			external_function_param_casadi_set_n_out(ext_fun_param_ptr, &sim_model_dyn_gnsf_get_matrices_fun_n_out);
			external_function_param_casadi_create(ext_fun_param_ptr, np);
			sim_in_set(config, dims, in, "gnsf_get_matrices_fun", ext_fun_param_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_param_ptr;
			mxSetField(plhs[0], 0, "dyn_gnsf_get_matrices_fun", tmp_mat);

			// gnsf_phi_fun
			ext_fun_param_ptr = (external_function_param_casadi *) malloc(1*sizeof(external_function_param_casadi));
			external_function_param_casadi_set_fun(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_fun);
			external_function_param_casadi_set_work(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_fun_work);
			external_function_param_casadi_set_sparsity_in(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_fun_sparsity_in);
			external_function_param_casadi_set_sparsity_out(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_fun_sparsity_out);
			external_function_param_casadi_set_n_in(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_fun_n_in);
			external_function_param_casadi_set_n_out(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_fun_n_out);
			external_function_param_casadi_create(ext_fun_param_ptr, np);
			sim_in_set(config, dims, in, "gnsf_phi_fun", ext_fun_param_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_param_ptr;
			mxSetField(plhs[0], 0, "dyn_gnsf_phi_fun", tmp_mat);

			// gnsf_phi_fun_jac_y
			ext_fun_param_ptr = (external_function_param_casadi *) malloc(1*sizeof(external_function_param_casadi));
			external_function_param_casadi_set_fun(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_fun_jac_y);
			external_function_param_casadi_set_work(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_fun_jac_y_work);
			external_function_param_casadi_set_sparsity_in(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_fun_jac_y_sparsity_in);
			external_function_param_casadi_set_sparsity_out(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_fun_jac_y_sparsity_out);
			external_function_param_casadi_set_n_in(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_fun_jac_y_n_in);
			external_function_param_casadi_set_n_out(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_fun_jac_y_n_out);
			external_function_param_casadi_create(ext_fun_param_ptr, np);
			sim_in_set(config, dims, in, "gnsf_phi_fun_jac_y", ext_fun_param_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_param_ptr;
			mxSetField(plhs[0], 0, "dyn_gnsf_phi_fun_jac_y", tmp_mat);

			// gnsf_phi_jac_y_uhat
			ext_fun_param_ptr = (external_function_param_casadi *) malloc(1*sizeof(external_function_param_casadi));
			external_function_param_casadi_set_fun(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_jac_y_uhat);
			external_function_param_casadi_set_work(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_jac_y_uhat_work);
			external_function_param_casadi_set_sparsity_in(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_jac_y_uhat_sparsity_in);
			external_function_param_casadi_set_sparsity_out(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_jac_y_uhat_sparsity_out);
			external_function_param_casadi_set_n_in(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_jac_y_uhat_n_in);
			external_function_param_casadi_set_n_out(ext_fun_param_ptr, &sim_model_dyn_gnsf_phi_jac_y_uhat_n_out);
			external_function_param_casadi_create(ext_fun_param_ptr, np);
			sim_in_set(config, dims, in, "gnsf_phi_jac_y_uhat", ext_fun_param_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_param_ptr;
			mxSetField(plhs[0], 0, "dyn_gnsf_phi_jac_y_uhat", tmp_mat);

			}
		else
			{
			// gnsf_f_lo_fun_jac_x1k1uz
			ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
			external_function_casadi_set_fun(ext_fun_ptr, &sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz);
			external_function_casadi_set_work(ext_fun_ptr, &sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz_work);
			external_function_casadi_set_sparsity_in(ext_fun_ptr, &sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz_sparsity_in);
			external_function_casadi_set_sparsity_out(ext_fun_ptr, &sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz_sparsity_out);
			external_function_casadi_set_n_in(ext_fun_ptr, &sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz_n_in);
			external_function_casadi_set_n_out(ext_fun_ptr, &sim_model_dyn_gnsf_f_lo_fun_jac_x1k1uz_n_out);
			external_function_casadi_create(ext_fun_ptr);
			sim_in_set(config, dims, in, "gnsf_f_lo_fun_jac_x1k1uz", ext_fun_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_ptr;
			mxSetField(plhs[0], 0, "dyn_gnsf_f_lo_fun_jac_x1k1uz", tmp_mat);

			// gnsf_get_matrices_fun
			ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
			external_function_casadi_set_fun(ext_fun_ptr, &sim_model_dyn_gnsf_get_matrices_fun);
			external_function_casadi_set_work(ext_fun_ptr, &sim_model_dyn_gnsf_get_matrices_fun_work);
			external_function_casadi_set_sparsity_in(ext_fun_ptr, &sim_model_dyn_gnsf_get_matrices_fun_sparsity_in);
			external_function_casadi_set_sparsity_out(ext_fun_ptr, &sim_model_dyn_gnsf_get_matrices_fun_sparsity_out);
			external_function_casadi_set_n_in(ext_fun_ptr, &sim_model_dyn_gnsf_get_matrices_fun_n_in);
			external_function_casadi_set_n_out(ext_fun_ptr, &sim_model_dyn_gnsf_get_matrices_fun_n_out);
			external_function_casadi_create(ext_fun_ptr);
			sim_in_set(config, dims, in, "gnsf_get_matrices_fun", ext_fun_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_ptr;
			mxSetField(plhs[0], 0, "dyn_gnsf_get_matrices_fun", tmp_mat);

			// gnsf_phi_fun
			ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
			external_function_casadi_set_fun(ext_fun_ptr, &sim_model_dyn_gnsf_phi_fun);
			external_function_casadi_set_work(ext_fun_ptr, &sim_model_dyn_gnsf_phi_fun_work);
			external_function_casadi_set_sparsity_in(ext_fun_ptr, &sim_model_dyn_gnsf_phi_fun_sparsity_in);
			external_function_casadi_set_sparsity_out(ext_fun_ptr, &sim_model_dyn_gnsf_phi_fun_sparsity_out);
			external_function_casadi_set_n_in(ext_fun_ptr, &sim_model_dyn_gnsf_phi_fun_n_in);
			external_function_casadi_set_n_out(ext_fun_ptr, &sim_model_dyn_gnsf_phi_fun_n_out);
			external_function_casadi_create(ext_fun_ptr);
			sim_in_set(config, dims, in, "gnsf_phi_fun", ext_fun_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_ptr;
			mxSetField(plhs[0], 0, "dyn_gnsf_phi_fun", tmp_mat);

			// gnsf_phi_fun_jac_y
			ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
			external_function_casadi_set_fun(ext_fun_ptr, &sim_model_dyn_gnsf_phi_fun_jac_y);
			external_function_casadi_set_work(ext_fun_ptr, &sim_model_dyn_gnsf_phi_fun_jac_y_work);
			external_function_casadi_set_sparsity_in(ext_fun_ptr, &sim_model_dyn_gnsf_phi_fun_jac_y_sparsity_in);
			external_function_casadi_set_sparsity_out(ext_fun_ptr, &sim_model_dyn_gnsf_phi_fun_jac_y_sparsity_out);
			external_function_casadi_set_n_in(ext_fun_ptr, &sim_model_dyn_gnsf_phi_fun_jac_y_n_in);
			external_function_casadi_set_n_out(ext_fun_ptr, &sim_model_dyn_gnsf_phi_fun_jac_y_n_out);
			external_function_casadi_create(ext_fun_ptr);
			sim_in_set(config, dims, in, "gnsf_phi_fun_jac_y", ext_fun_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_ptr;
			mxSetField(plhs[0], 0, "dyn_gnsf_phi_fun_jac_y", tmp_mat);

			// gnsf_phi_jac_y_uhat
			ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
			external_function_casadi_set_fun(ext_fun_ptr, &sim_model_dyn_gnsf_phi_jac_y_uhat);
			external_function_casadi_set_work(ext_fun_ptr, &sim_model_dyn_gnsf_phi_jac_y_uhat_work);
			external_function_casadi_set_sparsity_in(ext_fun_ptr, &sim_model_dyn_gnsf_phi_jac_y_uhat_sparsity_in);
			external_function_casadi_set_sparsity_out(ext_fun_ptr, &sim_model_dyn_gnsf_phi_jac_y_uhat_sparsity_out);
			external_function_casadi_set_n_in(ext_fun_ptr, &sim_model_dyn_gnsf_phi_jac_y_uhat_n_in);
			external_function_casadi_set_n_out(ext_fun_ptr, &sim_model_dyn_gnsf_phi_jac_y_uhat_n_out);
			external_function_casadi_create(ext_fun_ptr);
			sim_in_set(config, dims, in, "gnsf_phi_jac_y_uhat", ext_fun_ptr);
			// populate output struct
			tmp_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
			ptr = mxGetData(tmp_mat);
			ptr[0] = (long long) ext_fun_ptr;
			mxSetField(plhs[0], 0, "dyn_gnsf_phi_jac_y_uhat", tmp_mat);

			}
		}
	else
		{
		mexPrintf("\nsim_set_ext_fun_dyn_gnsf: method not supported %s\n", method);
		return;
		}
	
	return;

	}


