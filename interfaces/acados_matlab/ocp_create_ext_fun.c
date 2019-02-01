// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
//#include "acados/utils/external_function_generic.h"
//#include "acados_c/external_function_interface.h"
// mex
#include "mex.h"



// casadi functions for the model
//#include "model.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin sim_impl_ext_fun_create\n");

	/* RHS */



	/* LHS */

	// field names of output struct
	char *fieldnames[12];
	fieldnames[0] = (char*)mxMalloc(50);
	fieldnames[1] = (char*)mxMalloc(50);
	fieldnames[2] = (char*)mxMalloc(50);
	fieldnames[3] = (char*)mxMalloc(50);
	fieldnames[4] = (char*)mxMalloc(50);
	fieldnames[5] = (char*)mxMalloc(50);
	fieldnames[6] = (char*)mxMalloc(50);
	fieldnames[7] = (char*)mxMalloc(50);
	fieldnames[8] = (char*)mxMalloc(50);
	fieldnames[9] = (char*)mxMalloc(50);
	fieldnames[10] = (char*)mxMalloc(50);
	fieldnames[11] = (char*)mxMalloc(50);

	memcpy(fieldnames[0],"expl_ode_fun",sizeof("expl_ode_fun"));
	memcpy(fieldnames[1],"expl_vde_for",sizeof("expl_vde_for"));
	memcpy(fieldnames[2],"expl_vde_adj",sizeof("expl_vde_adj"));
	memcpy(fieldnames[3],"impl_ode_fun",sizeof("impl_ode_fun"));
	memcpy(fieldnames[4],"impl_ode_fun_jac_x_xdot",sizeof("impl_ode_fun_jac_x_xdot"));
	memcpy(fieldnames[5],"impl_ode_jac_x_xdot_u",sizeof("impl_ode_jac_x_xdot_u"));
	memcpy(fieldnames[6],"h_fun_jac_ut_xt",sizeof("h_fun_jac_ut_xt"));
	memcpy(fieldnames[7],"h_e_fun_jac_ut_xt",sizeof("h_e_fun_jac_ut_xt"));
	memcpy(fieldnames[8],"y_fun_jac_ut_xt",sizeof("y_fun_jac_ut_xt"));
	memcpy(fieldnames[9],"y_e_fun_jac_ut_xt",sizeof("y_e_fun_jac_ut_xt"));
	memcpy(fieldnames[10],"ext_cost_jac_hes",sizeof("ext_cost_jac_hes"));
	memcpy(fieldnames[11],"ext_cost_e_jac_hes",sizeof("ext_cost_e_jac_hes"));

	// create output struct
	plhs[0] = mxCreateStructMatrix(1, 1, 12, (const char **) fieldnames);

	mxFree( fieldnames[0] );
	mxFree( fieldnames[1] );
	mxFree( fieldnames[2] );
	mxFree( fieldnames[3] );
	mxFree( fieldnames[4] );
	mxFree( fieldnames[5] );
	mxFree( fieldnames[6] );
	mxFree( fieldnames[7] );
	mxFree( fieldnames[8] );
	mxFree( fieldnames[9] );
	mxFree( fieldnames[10] );
	mxFree( fieldnames[11] );



	return;

	}




