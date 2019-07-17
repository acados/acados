// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados/utils/external_function_generic.h"
#include "acados_c/external_function_interface.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin sim_ext_fun_destroy\n");

	int ii;
	long long *ptr;



	/* RHS */

	// model

	// C_sim_ext_fun

	external_function_param_casadi *ext_fun_param_ptr;
	int struct_size = mxGetNumberOfFields( prhs[1] );
	for(ii=0; ii<struct_size; ii++)
		{
//		printf("\n%s\n", mxGetFieldNameByNumber( prhs[1], ii) );
		ptr = (long long *) mxGetData( mxGetFieldByNumber( prhs[1], 0, ii ) );
		ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
		if(ext_fun_param_ptr!=0)
			{
			external_function_param_casadi_free(ext_fun_param_ptr);
			free(ext_fun_param_ptr);
			}
		}

	/* return */

	return;

	}

