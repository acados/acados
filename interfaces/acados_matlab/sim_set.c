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



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin sim_set\n");

	int ii;
	long long *ptr;

	/* RHS */

	// model

	// opts
	char *method = mxArrayToString( mxGetField( prhs[1], 0, "method" ) );

	// C_sim

	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "config" ) );
	sim_config *config = (sim_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dims" ) );
	void *dims = (void *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "in" ) );
	sim_in *in = (sim_in *) ptr[0];

	// field
	char *field = mxArrayToString( prhs[4] );
//	mexPrintf("\n%s\n", field);


	// value
	if(!strcmp(field, "T"))
		{
		double *T = mxGetPr( prhs[5] );
		sim_in_set(config, dims, in, "T", T);
		}
	else if(!strcmp(field, "x"))
		{
		double *x = mxGetPr( prhs[5] );
		sim_in_set(config, dims, in, "x", x);
		}
	else if(!strcmp(field, "u"))
		{
		double *u = mxGetPr( prhs[5] );
		sim_in_set(config, dims, in, "u", u);
		}
	else if(!strcmp(field, "p"))
		{
		double *p = mxGetPr( prhs[5] );
		external_function_param_casadi *ext_fun_param_ptr;
		int struct_size = mxGetNumberOfFields( prhs[3] );
		for(ii=0; ii<struct_size; ii++)
			{
//			printf("\n%s\n", mxGetFieldNameByNumber( prhs[3], ii) );
			ptr = (long long *) mxGetData( mxGetFieldByNumber( prhs[3], 0, ii ) );
			if(ptr[0]!=0)
				{
				ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
				ext_fun_param_ptr->set_param(ext_fun_param_ptr, p);
				}
			}
		}
	else if(!strcmp(field, "seed_adj"))
		{
		double *seed_adj = mxGetPr( prhs[5] );
		sim_in_set(config, dims, in, "seed_adj", seed_adj);
		}
	else
		{
		mexPrintf("\nsim_set: field not supported: %s\n", field);
		return;
		}



	/* return */

	return;

	}


