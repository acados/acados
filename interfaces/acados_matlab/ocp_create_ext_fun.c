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
	char *fieldnames[18];
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
	fieldnames[12] = (char*)mxMalloc(50);
	fieldnames[13] = (char*)mxMalloc(50);
	fieldnames[14] = (char*)mxMalloc(50);
	fieldnames[15] = (char*)mxMalloc(50);
	fieldnames[16] = (char*)mxMalloc(50);
	fieldnames[17] = (char*)mxMalloc(50);

	memcpy(fieldnames[0],"dyn_expl_ode_fun",sizeof("dyn_expl_ode_fun"));
	memcpy(fieldnames[1],"dyn_expl_vde_for",sizeof("dyn_expl_vde_for"));
	memcpy(fieldnames[2],"dyn_expl_vde_adj",sizeof("dyn_expl_vde_adj"));
	memcpy(fieldnames[3],"dyn_expl_ode_hes",sizeof("dyn_expl_ode_hes"));
	memcpy(fieldnames[4],"dyn_impl_ode_fun",sizeof("dyn_impl_ode_fun"));
	memcpy(fieldnames[5],"dyn_impl_ode_fun_jac_x_xdot",sizeof("dyn_impl_ode_fun_jac_x_xdot"));
	memcpy(fieldnames[6],"dyn_impl_ode_jac_x_xdot_u",sizeof("dyn_impl_ode_jac_x_xdot_u"));
	memcpy(fieldnames[7],"dyn_impl_ode_hess",sizeof("dyn_impl_ode_hess"));
	memcpy(fieldnames[8],"dyn_disc_phi_fun_jac",sizeof("dyn_disc_phi_fun_jac"));
	memcpy(fieldnames[9],"dyn_disc_phi_fun_jac_hess",sizeof("dyn_disc_phi_fun_jac_hess"));
	memcpy(fieldnames[10],"constr_h_fun_jac_ut_xt",sizeof("constr_h_fun_jac_ut_xt"));
	memcpy(fieldnames[11],"constr_h_fun_jac_ut_xt_hess",sizeof("constr_h_fun_jac_ut_xt_hess"));
	memcpy(fieldnames[12],"constr_h_e_fun_jac_ut_xt",sizeof("constr_h_e_fun_jac_ut_xt"));
	memcpy(fieldnames[13],"constr_h_e_fun_jac_ut_xt_hess",sizeof("constr_h_e_fun_jac_ut_xt_hess"));
	memcpy(fieldnames[14],"cost_y_fun_jac_ut_xt",sizeof("cost_y_fun_jac_ut_xt"));
	memcpy(fieldnames[15],"cost_y_e_fun_jac_ut_xt",sizeof("cost_y_e_fun_jac_ut_xt"));
	memcpy(fieldnames[16],"cost_ext_cost_jac_hes",sizeof("cost_ext_cost_jac_hes"));
	memcpy(fieldnames[17],"cost_ext_cost_e_jac_hes",sizeof("cost_ext_cost_e_jac_hes"));

	// create output struct
	plhs[0] = mxCreateStructMatrix(1, 1, 18, (const char **) fieldnames);

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
	mxFree( fieldnames[12] );
	mxFree( fieldnames[13] );
	mxFree( fieldnames[14] );
	mxFree( fieldnames[15] );
	mxFree( fieldnames[16] );
	mxFree( fieldnames[17] );



	return;

	}




