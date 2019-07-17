// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados_c/ocp_nlp_interface.h"
#include "acados/utils/external_function_generic.h"
#include "acados_c/external_function_interface.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin ocp_ext_fun_destroy\n");

	int ii, jj, kk;
	long long *ptr;
	int Nf;
	mxArray *mex_field;

	/* RHS */

	// model

	// C_ocp

	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dims" ) );
	ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];

	int N = dims->N;

	// XXX hard-code number and size of phases for now !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	int NN[] = {N, 1};


	// C_ocp_ext_fun

	//
	external_function_param_casadi *ext_fun_param_ptr;
	int struct_size = mxGetNumberOfFields( prhs[2] );
	for(ii=0; ii<struct_size; ii++)
		{
//		printf("\n%s\n", mxGetFieldNameByNumber( prhs[2], ii) );
		mex_field = mxGetFieldByNumber( prhs[2], 0, ii );
		ptr = (long long *) mxGetData( mex_field );
		Nf = mxGetN( mex_field );
		for(jj=0; jj<Nf; jj++)
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[jj];
			if(ext_fun_param_ptr!=0)
				{
				for(kk=0; kk<NN[jj]; kk++)
					{
					external_function_param_casadi_free(ext_fun_param_ptr+kk);
					}
				free(ext_fun_param_ptr);
				}
			}
		}

	/* return */

	return;

	}


