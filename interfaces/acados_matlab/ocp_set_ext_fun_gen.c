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


// macro to string
#define STR(x) STR_AGAIN(x)
#define STR_AGAIN(x) #x

// glue macros
#define GLUE2(x,y) GLUE2_AGAIN(x,y)
#define GLUE2_AGAIN(x,y) x##y

// macro bricks
#define WORK _work
#define SP_IN _sparsity_in
#define SP_OUT _sparsity_out
#define N_IN _n_in
#define N_OUT _n_out


// casadi functions for the model
int FUN_NAME(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int GLUE2(FUN_NAME,WORK)(int *, int *, int *, int *);
const int *GLUE2(FUN_NAME,SP_IN)(int);
const int *GLUE2(FUN_NAME,SP_OUT)(int);
int GLUE2(FUN_NAME,N_IN)();
int GLUE2(FUN_NAME,N_OUT)();



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin ocp_expl_ext_fun_create\n");

	// sizeof(long long) == sizeof(void *) = 64 !!!
	int ii, status;
	long long *ptr;
//	mxInt64 *ptr;



	/* RHS */

	// C_ocp

	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "config" ) );
	ocp_nlp_config *config = (ocp_nlp_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "dims" ) );
	ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "in" ) );
	ocp_nlp_in *in = (ocp_nlp_in *) ptr[0];

	// C_ocp_ext_fun

	// model

	int np = 0;
	if(mxGetField( prhs[2], 0, "dim_np" )!=NULL) // TODO bool
		{
		np = mxGetScalar( mxGetField( prhs[2], 0, "dim_np" ) );
		}


	
	/* LHS */

	/* copy existing fields */

//	plhs[0] = mxCreateSharedDataCopy(prhs[1]);
	plhs[0] = mxDuplicateArray(prhs[1]);



	/* populate new fields */

	external_function_casadi *ext_fun_ptr;
	external_function_param_casadi *ext_fun_param_ptr;

	mxArray *tmp_mat;

	if(!strcmp(STR(PARAM), "true")) // TODO bool
		{
		// expl_ode_fun
		ext_fun_param_ptr = (external_function_param_casadi *) malloc((N1-N0+1)*sizeof(external_function_param_casadi));
		for(ii=0; ii<N1-N0+1; ii++)
			{
			external_function_param_casadi_set_fun(ext_fun_param_ptr+ii, &FUN_NAME);
			external_function_param_casadi_set_work(ext_fun_param_ptr+ii, &GLUE2(FUN_NAME,WORK));
			external_function_param_casadi_set_sparsity_in(ext_fun_param_ptr+ii, &GLUE2(FUN_NAME,SP_IN));
			external_function_param_casadi_set_sparsity_out(ext_fun_param_ptr+ii, &GLUE2(FUN_NAME,SP_OUT));
			external_function_param_casadi_set_n_in(ext_fun_param_ptr+ii, &GLUE2(FUN_NAME,N_IN));
			external_function_param_casadi_set_n_out(ext_fun_param_ptr+ii, &GLUE2(FUN_NAME,N_OUT));
			external_function_param_casadi_create(ext_fun_param_ptr+ii, np);
			status = SETTER(config, dims, in, N0+ii, STR(SET_FIELD), ext_fun_param_ptr+ii);
			}
		// populate output struct
		ptr = mxGetData(mxGetField(plhs[0], 0, STR(MEX_FIELD)));
		ptr[PHASE] = (long long) ext_fun_param_ptr;
		}
	else
		{
		// expl_ode_fun
		ext_fun_ptr = (external_function_casadi *) malloc((N1-N0+1)*sizeof(external_function_casadi));
		for(ii=0; ii<N1-N0+1; ii++)
			{
			external_function_casadi_set_fun(ext_fun_ptr+ii, &FUN_NAME);
			external_function_casadi_set_work(ext_fun_ptr+ii, &GLUE2(FUN_NAME,WORK));
			external_function_casadi_set_sparsity_in(ext_fun_ptr+ii, &GLUE2(FUN_NAME,SP_IN));
			external_function_casadi_set_sparsity_out(ext_fun_ptr+ii, &GLUE2(FUN_NAME,SP_OUT));
			external_function_casadi_set_n_in(ext_fun_ptr+ii, &GLUE2(FUN_NAME,N_IN));
			external_function_casadi_set_n_out(ext_fun_ptr+ii, &GLUE2(FUN_NAME,N_OUT));
			external_function_casadi_create(ext_fun_ptr+ii);
			status = SETTER(config, dims, in, N0+ii, STR(SET_FIELD), ext_fun_ptr+ii);
			}
		// populate output struct
		ptr = mxGetData(mxGetField(plhs[0], 0, STR(MEX_FIELD)));
		ptr[PHASE] = (long long) ext_fun_ptr;
		}
	
	return;

	}





