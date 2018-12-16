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

//	mexPrintf("\nin ocp_expl_ext_fun_create\n");

	// sizeof(long long) == sizeof(void *) = 64 !!!
	long long *ptr;



	/* RHS */

	// C_ocp_ext_fun

	// opts

	// TODO use them !!!
//	bool sens_forw = mxGetScalar( mxGetField( prhs[0], 0, "sens_forw" ) );
//	mexPrintf("\n%d\n", sens_forw);
	char *sim_solver = mxArrayToString( mxGetField( prhs[1], 0, "sim_solver" ) );
//	mexPrintf("\n%s\n", sim_solver);



	/* LHS */





	/* populate input struc */

	external_function_casadi *ext_fun_ptr;

	// TODO templetize the casadi function names !!!
	if(!strcmp(sim_solver, "erk"))
		{
		// expl_ode_fun
		ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
		ext_fun_ptr->casadi_fun = &model_expl_ode_fun;
		ext_fun_ptr->casadi_work = &model_expl_ode_fun_work;
		ext_fun_ptr->casadi_sparsity_in = &model_expl_ode_fun_sparsity_in;
		ext_fun_ptr->casadi_sparsity_out = &model_expl_ode_fun_sparsity_out;
		ext_fun_ptr->casadi_n_in = &model_expl_ode_fun_n_in;
		ext_fun_ptr->casadi_n_out = &model_expl_ode_fun_n_out;
		external_function_casadi_create(ext_fun_ptr);
		// populate output struct
		mxArray *expl_ode_fun_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
		ptr = mxGetData(expl_ode_fun_mat);
		ptr[0] = (long long) ext_fun_ptr;
		mxSetField((mxArray*) prhs[0], 0, "expl_ode_fun", expl_ode_fun_mat);

		// expl_vde_for
		ext_fun_ptr = (external_function_casadi *) malloc(1*sizeof(external_function_casadi));
		ext_fun_ptr->casadi_fun = &model_expl_vde_for;
		ext_fun_ptr->casadi_work = &model_expl_vde_for_work;
		ext_fun_ptr->casadi_sparsity_in = &model_expl_vde_for_sparsity_in;
		ext_fun_ptr->casadi_sparsity_out = &model_expl_vde_for_sparsity_out;
		ext_fun_ptr->casadi_n_in = &model_expl_vde_for_n_in;
		ext_fun_ptr->casadi_n_out = &model_expl_vde_for_n_out;
		external_function_casadi_create(ext_fun_ptr);
		// populate output struct
		mxArray *expl_vde_for_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
		ptr = mxGetData(expl_vde_for_mat);
		ptr[0] = (long long) ext_fun_ptr;
		mxSetField((mxArray*) prhs[0], 0, "expl_vde_for", expl_vde_for_mat);
		}
	else
		{
		mexPrintf("\nocp_set_ext_fun_expl: sim_solver not supported %s\n", sim_solver);
		return;
		}
	
	return;

	}




