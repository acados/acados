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
#include "model.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin sim_impl_ext_fun_create\n");

	// sizeof(long long) == sizeof(void *) = 64 !!!
	long long *ptr;



	/* RHS */

	// opts

//	bool sens_forw = mxGetScalar( mxGetField( prhs[0], 0, "sens_forw" ) );
//	mexPrintf("\n%d\n", sens_forw);
//	char *sim_solver = mxArrayToString( mxGetField( prhs[0], 0, "sim_solver" ) );
//	mexPrintf("\n%s\n", sim_solver);



	/* LHS */

#if 0
	// field names of output struct
	char *fieldnames[7];
	fieldnames[0] = (char*)mxMalloc(50);
	fieldnames[1] = (char*)mxMalloc(50);
	fieldnames[2] = (char*)mxMalloc(50);
	fieldnames[3] = (char*)mxMalloc(50);
	fieldnames[4] = (char*)mxMalloc(50);
	fieldnames[5] = (char*)mxMalloc(50);
	fieldnames[6] = (char*)mxMalloc(50);

	memcpy(fieldnames[0],"expl_ode_fun",sizeof("expl_ode_fun"));
	memcpy(fieldnames[1],"expl_vde_for",sizeof("expl_vde_for"));
	memcpy(fieldnames[2],"expl_vde_adj",sizeof("expl_vde_adj"));
	memcpy(fieldnames[3],"impl_ode_fun",sizeof("impl_ode_fun"));
	memcpy(fieldnames[4],"impl_ode_fun_jac_x_xdot",sizeof("impl_ode_fun_jac_x_xdot"));
	memcpy(fieldnames[5],"impl_ode_jac_x_xdot_u",sizeof("impl_ode_jac_x_xdot_u"));
	memcpy(fieldnames[6],"h_fun_jac_ut_xt",sizeof("h_fun_jac_ut_xt"));

	// create output struct
	plhs[0] = mxCreateStructMatrix(1, 1, 7, (const char **) fieldnames);

	mxFree( fieldnames[0] );
	mxFree( fieldnames[1] );
	mxFree( fieldnames[2] );
	mxFree( fieldnames[3] );
	mxFree( fieldnames[4] );
	mxFree( fieldnames[5] );
	mxFree( fieldnames[6] );
#endif

	

	
	// TODO do it properly with nrhs
	// XXX assume the struct is already created, and just passed over from RHS to LHS
	plhs[0] = prhs[0];



	external_function_casadi *ext_fun_ptr;



	// TODO templetize the casadi function names !!!
//	if(!strcmp(sim_solver, "irk"))
//		{
		// h_fun_jac_ut_xt
		ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
		ext_fun_ptr->casadi_fun = &model_h_fun_jac_ut_xt;
		ext_fun_ptr->casadi_work = &model_h_fun_jac_ut_xt_work;
		ext_fun_ptr->casadi_sparsity_in = &model_h_fun_jac_ut_xt_sparsity_in;
		ext_fun_ptr->casadi_sparsity_out = &model_h_fun_jac_ut_xt_sparsity_out;
		ext_fun_ptr->casadi_n_in = &model_h_fun_jac_ut_xt_n_in;
		ext_fun_ptr->casadi_n_out = &model_h_fun_jac_ut_xt_n_out;
		external_function_casadi_create(ext_fun_ptr);
		// populate output struct
		mxArray *h_fun_jac_ut_xt_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
		ptr = mxGetData(h_fun_jac_ut_xt_mat);
		ptr[0] = (long long) ext_fun_ptr;
		mxSetField(plhs[0], 0, "h_fun_jac_ut_xt", h_fun_jac_ut_xt_mat);
//		}
//	else
//		{
//		mexPrintf("\nocp_impl_ext_fun_create: sim_solver not supported %s\n", sim_solver);
//		return;
//		}
	
	return;

	}




