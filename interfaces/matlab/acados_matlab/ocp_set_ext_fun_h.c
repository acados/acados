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
#include "ocp_model.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin sim_impl_ext_fun_create\n");

	// sizeof(long long) == sizeof(void *) = 64 !!!
	long long *ptr;



	/* RHS */

	// C_ocp_ext_fun

	// model

	int nh;			bool set_h = false;

	if (mxGetField( prhs[1], 0, "nh" )!=NULL)
		{
		set_h = true;
		nh = mxGetScalar( mxGetField( prhs[1], 0, "nh" ) );
		}

	// opts

//	bool sens_forw = mxGetScalar( mxGetField( prhs[0], 0, "sens_forw" ) );
//	mexPrintf("\n%d\n", sens_forw);
//	char *sim_method = mxArrayToString( mxGetField( prhs[2], 0, "sim_method" ) );
//	mexPrintf("\n%s\n", sim_method);



	/* LHS */




	/* populate input struc */

	external_function_casadi *ext_fun_ptr;

	// TODO templetize the casadi function names !!!
	if(set_h && nh>0)
		{
		// h_fun_jac_ut_xt
		ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
		ext_fun_ptr->casadi_fun = &ocp_model_h_fun_jac_ut_xt;
		ext_fun_ptr->casadi_work = &ocp_model_h_fun_jac_ut_xt_work;
		ext_fun_ptr->casadi_sparsity_in = &ocp_model_h_fun_jac_ut_xt_sparsity_in;
		ext_fun_ptr->casadi_sparsity_out = &ocp_model_h_fun_jac_ut_xt_sparsity_out;
		ext_fun_ptr->casadi_n_in = &ocp_model_h_fun_jac_ut_xt_n_in;
		ext_fun_ptr->casadi_n_out = &ocp_model_h_fun_jac_ut_xt_n_out;
		external_function_casadi_create(ext_fun_ptr);
		// populate output struct
		mxArray *h_fun_jac_ut_xt_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
		ptr = mxGetData(h_fun_jac_ut_xt_mat);
		ptr[0] = (long long) ext_fun_ptr;
		mxSetField((mxArray*) prhs[0], 0, "h_fun_jac_ut_xt", h_fun_jac_ut_xt_mat);
		}
	
	return;

	}




