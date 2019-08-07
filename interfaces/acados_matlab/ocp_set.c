// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados_c/ocp_nlp_interface.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin ocp_get\n");

	long long *ptr;
	int ii, jj, kk;
	int Nf;
	mxArray *mex_field;

	/* RHS */

	// model

	// C_ocp

	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "config" ) );
	ocp_nlp_config *config = (ocp_nlp_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dims" ) );
	ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "in" ) );
	ocp_nlp_in *in = (ocp_nlp_in *) ptr[0];
	// out
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "out" ) );
	ocp_nlp_out *out = (ocp_nlp_out *) ptr[0];

	// opts
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "opts" ) );
	void *opts = (void *) ptr[0];

	// field
	char *field = mxArrayToString( prhs[4] );
//	mexPrintf("\n%s\n", field);



#if 0
	char *module = NULL;

	char *char_dot = strchr(field, '.');
	if(char_dot!=NULL)
		{
		int size = char_dot-field;
		mexPrintf("\nyep dot found %d\n", size);
		module = malloc((size+1)*sizeof(char));
//		strncpy(field, module, size);
		for(ii=0; ii<size; ii++)
			{
			module[ii] = field[ii];
			}
		module[size] = '\0'; // add end of string
		mexPrintf("\nmodule: %s\n", module);
		}
	else
		{
		mexPrintf("\nno dot found\n");
		}
#endif



	int N = dims->N;
	int nu = dims->nu[0];
	int nx = dims->nx[0];

	// XXX hard-code number and size of phases for now !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	int NN[] = {N, 1};


	// TODO implement with LHS ?????
	// value
	if (!strcmp(field, "constr_x0"))
		{
		double *x0 = mxGetPr( prhs[5] );
		ocp_nlp_constraints_model_set(config, dims, in, 0, "lbx", x0);
		ocp_nlp_constraints_model_set(config, dims, in, 0, "ubx", x0);
		}
	else if (!strcmp(field, "cost_yr"))
		{
		if(nrhs==6)
			{
			double *yr = mxGetPr( prhs[5] );
			for (ii=0; ii<N; ii++)
				{
				ocp_nlp_cost_model_set(config, dims, in, ii, "y_ref", yr);
				}
			}
		else if(nrhs==7)
			{
			double *yr = mxGetPr( prhs[5] );
			int stage = mxGetScalar( prhs[6] );
			ocp_nlp_cost_model_set(config, dims, in, stage, "y_ref", yr);
			}
		else
			{
			mexPrintf("\nocp_set: wrong nrhs: %d\n", nrhs);
			goto end;
			}
		}
	else if (!strcmp(field, "cost_yr_e"))
		{
		double *yr_e = mxGetPr( prhs[5] );
		ocp_nlp_cost_model_set(config, dims, in, N, "y_ref", yr_e);
		}
	else if (!strcmp(field, "init_x"))
		{
		double *init_x = mxGetPr( prhs[5] );
		for (ii=0; ii<=N; ii++)
			{
			ocp_nlp_out_set(config, dims, out, ii, "x", init_x+ii*nx);
			}
		}
	else if (!strcmp(field, "init_u"))
		{
		double *init_u = mxGetPr( prhs[5] );
		for (ii=0; ii<N; ii++)
			{
			ocp_nlp_out_set(config, dims, out, ii, "u", init_u+ii*nu);
			}
		}
	else if (!strcmp(field, "init_pi"))
		{
		double *init_pi = mxGetPr( prhs[5] );
		for (ii=0; ii<N; ii++)
			{
			ocp_nlp_out_set(config, dims, out, ii, "pi", init_pi+ii*nx);
			}
		}
	// TODO make setters for ALL numerical data
	else if(!strcmp(field, "p"))
		{
		double *p = mxGetPr( prhs[5] );
		external_function_param_casadi *ext_fun_param_ptr;
		int struct_size = mxGetNumberOfFields( prhs[3] );
		for(ii=0; ii<struct_size; ii++)
			{
//			printf("\n%s\n", mxGetFieldNameByNumber( prhs[3], ii) );
			mex_field = mxGetFieldByNumber( prhs[3], 0, ii );
			ptr = (long long *) mxGetData( mex_field );
			Nf = mxGetN( mex_field );
			int Nf_sum = 0;
			for(jj=0; jj<Nf; jj++)
				{
				ext_fun_param_ptr = (external_function_param_casadi *) ptr[jj];
				if(ext_fun_param_ptr!=0)
					{
					if(nrhs==6)
						{
						for(kk=0; kk<NN[jj]; kk++)
							{
							(ext_fun_param_ptr+kk)->set_param(ext_fun_param_ptr+kk, p);
							}
						}
					else if(nrhs==7)
						{
						int stage = mxGetScalar( prhs[6] );
						if(stage>=Nf_sum & stage<Nf_sum+NN[jj])
							{
							(ext_fun_param_ptr+stage)->set_param(ext_fun_param_ptr+stage, p);
							}
						}
					else
						{
						mexPrintf("\nocp_set: wrong nrhs: %d\n", nrhs);
						goto end;
						}
					}
				Nf_sum += NN[jj];
				}
			}
		}
	else if (!strcmp(field, "nlp_solver_max_iter"))
	{
		double *max_it_double = mxGetPr( prhs[5] );
		int nlp_solver_max_iter = max_it_double[0];
		ocp_nlp_opts_set(config, opts, "max_iter", &nlp_solver_max_iter);
	}
	else
		{
		mexPrintf("\nocp_set: field not supported: %s\n", field);
		goto end;
		}
	

	/* return */
end:
#if 0
	if(module!=NULL)
		{
		free(module);
		}
#endif

	return;

	}




