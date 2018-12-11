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

//	mexPrintf("\nin sim_ext_fun_create\n");

	// sizeof(long long) == sizeof(void *) = 64 !!!
	long long *ptr;



	/* RHS */

	// opts

	bool sens_forw = mxGetScalar( mxGetField( prhs[0], 0, "sens_forw" ) );
//	mexPrintf("\n%d\n", sens_forw);
	char *scheme = mxArrayToString( mxGetField( prhs[0], 0, "scheme" ) );
//	mexPrintf("\n%s\n", scheme);



	/* LHS */

	// field names of output struct
	char *fieldnames[6];
	fieldnames[0] = (char*)mxMalloc(50);
	fieldnames[1] = (char*)mxMalloc(50);
	fieldnames[2] = (char*)mxMalloc(50);
	fieldnames[3] = (char*)mxMalloc(50);
	fieldnames[4] = (char*)mxMalloc(50);
	fieldnames[5] = (char*)mxMalloc(50);

	memcpy(fieldnames[0],"expl_ode_fun",sizeof("expl_ode_fun"));
	memcpy(fieldnames[1],"expl_vde_for",sizeof("expl_vde_for"));
	memcpy(fieldnames[2],"expl_vde_adj",sizeof("expl_vde_adj"));
	memcpy(fieldnames[3],"impl_ode_fun",sizeof("impl_ode_fun"));
	memcpy(fieldnames[4],"impl_ode_fun_jac_x_xdot",sizeof("impl_ode_fun_jac_x_xdot"));
	memcpy(fieldnames[5],"impl_ode_jac_x_xdot_u",sizeof("impl_ode_jac_x_xdot_u"));

	// create output struct
	plhs[0] = mxCreateStructMatrix(1, 1, 6, (const char **) fieldnames);

	mxFree( fieldnames[0] );
	mxFree( fieldnames[1] );
	mxFree( fieldnames[2] );
	mxFree( fieldnames[3] );
	mxFree( fieldnames[4] );
	mxFree( fieldnames[5] );



	external_function_casadi *expl_ode_fun,
	                         *expl_vde_for,
	                         *expl_vde_adj,
	                         *impl_ode_fun,
	                         *impl_ode_fun_jac_x_xdot,
	                         *impl_ode_jac_x_xdot_u;



	// TODO templetize the casadi function names !!!
	if(!strcmp(scheme, "irk"))
		{
		// impl_ode_fun
		impl_ode_fun = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
		impl_ode_fun->casadi_fun = &model_impl_ode_fun;
		impl_ode_fun->casadi_work = &model_impl_ode_fun_work;
		impl_ode_fun->casadi_sparsity_in = &model_impl_ode_fun_sparsity_in;
		impl_ode_fun->casadi_sparsity_out = &model_impl_ode_fun_sparsity_out;
		impl_ode_fun->casadi_n_in = &model_impl_ode_fun_n_in;
		impl_ode_fun->casadi_n_out = &model_impl_ode_fun_n_out;
		external_function_casadi_create(impl_ode_fun);
		// populate output struct
		mxArray *impl_ode_fun_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
		ptr = mxGetData(impl_ode_fun_mat);
		ptr[0] = (long long) impl_ode_fun;
		mxSetField(plhs[0], 0, "impl_ode_fun", impl_ode_fun_mat);

		// impl_ode_fun_jac_x_xdot
		impl_ode_fun_jac_x_xdot = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
		impl_ode_fun_jac_x_xdot->casadi_fun = &model_impl_ode_fun_jac_x_xdot;
		impl_ode_fun_jac_x_xdot->casadi_work = &model_impl_ode_fun_jac_x_xdot_work;
		impl_ode_fun_jac_x_xdot->casadi_sparsity_in = &model_impl_ode_fun_jac_x_xdot_sparsity_in;
		impl_ode_fun_jac_x_xdot->casadi_sparsity_out = &model_impl_ode_fun_jac_x_xdot_sparsity_out;
		impl_ode_fun_jac_x_xdot->casadi_n_in = &model_impl_ode_fun_jac_x_xdot_n_in;
		impl_ode_fun_jac_x_xdot->casadi_n_out = &model_impl_ode_fun_jac_x_xdot_n_out;
		external_function_casadi_create(impl_ode_fun_jac_x_xdot);
		// populate output struct
		mxArray *impl_ode_fun_jac_x_xdot_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
		ptr = mxGetData(impl_ode_fun_jac_x_xdot_mat);
		ptr[0] = (long long) impl_ode_fun_jac_x_xdot;
		mxSetField(plhs[0], 0, "impl_ode_fun_jac_x_xdot", impl_ode_fun_jac_x_xdot_mat);

		// impl_ode_jac_x_xdot_u
		impl_ode_jac_x_xdot_u = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
		impl_ode_jac_x_xdot_u->casadi_fun = &model_impl_ode_jac_x_xdot_u;
		impl_ode_jac_x_xdot_u->casadi_work = &model_impl_ode_jac_x_xdot_u_work;
		impl_ode_jac_x_xdot_u->casadi_sparsity_in = &model_impl_ode_jac_x_xdot_u_sparsity_in;
		impl_ode_jac_x_xdot_u->casadi_sparsity_out = &model_impl_ode_jac_x_xdot_u_sparsity_out;
		impl_ode_jac_x_xdot_u->casadi_n_in = &model_impl_ode_jac_x_xdot_u_n_in;
		impl_ode_jac_x_xdot_u->casadi_n_out = &model_impl_ode_jac_x_xdot_u_n_out;
		external_function_casadi_create(impl_ode_jac_x_xdot_u);
		// populate output struct
		mxArray *impl_ode_jac_x_xdot_u_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
		ptr = mxGetData(impl_ode_jac_x_xdot_u_mat);
		ptr[0] = (long long) impl_ode_jac_x_xdot_u;
		mxSetField(plhs[0], 0, "impl_ode_jac_x_xdot_u", impl_ode_jac_x_xdot_u_mat);
		}
	else
		{
		mexPrintf("\nsim_impl_ext_fun_create: scheme not supported %s\n", scheme);
		return;
		}
	
	return;

	}


