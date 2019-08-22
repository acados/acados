/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */

// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados_c/ocp_nlp_interface.h"
// mex
#include "mex.h"
#include "mex_macros.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	long long *ptr;
	int ii, jj, kk, acados_size;
	mxArray *mex_field;
	char fun_name[50] = "ocp_set";
	char buffer [100]; // for error messages


	/* RHS */
	// config
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "config" ) );
	ocp_nlp_config *config = (ocp_nlp_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dims" ) );
	ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];
	// opts
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "opts" ) );
	void *opts = (void *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "in" ) );
	ocp_nlp_in *in = (ocp_nlp_in *) ptr[0];
	// out
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "out" ) );
	ocp_nlp_out *out = (ocp_nlp_out *) ptr[0];

	const mxArray *C_ext_fun_pointers = prhs[3];
	// field
	char *field = mxArrayToString( prhs[4] );
	// value
	double *value = mxGetPr( prhs[5] );
	int matlab_size = (int) mxGetNumberOfElements( prhs[5] );

	int N = dims->N;
	int nu = dims->nu[0];
	int nx = dims->nx[0];

	// TODO(oj): Can this be removed?
#if 0
	char *module = NULL;

	char *char_dot = strchr(field, '.');
	if (char_dot!=NULL)
		{
		int size = char_dot-field;
		mexPrintf("\nyep dot found %d\n", size);
		module = malloc((size+1)*sizeof(char));
//		strncpy(field, module, size);
		for (ii=0; ii<size; ii++)
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

	// stage
	int s0, se;
	if (nrhs==6)
	{
		s0 = 0;
		se = N;
	}
	else if (nrhs==7)
	{
		s0 = mxGetScalar( prhs[6] );
		se = s0 + 1;
	}
	else
	{
		sprintf(buffer, "ocp_set: wrong nrhs: %d\n", nrhs);
		mexErrMsgTxt(buffer);
	}

	// XXX hard-code number and size of phases for now !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	int NN[] = {N, 1}; // number of phases, i.e. shooting nodes with same dimensions
	int Nf = 2; // number of phases
	// Nf = mxGetN( mex_field );

	// value
	if (!strcmp(field, "constr_x0"))
	{
		acados_size = nx;
		MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);
		ocp_nlp_constraints_model_set(config, dims, in, 0, "lbx", value);
		ocp_nlp_constraints_model_set(config, dims, in, 0, "ubx", value);
	}
	else if (!strcmp(field, "cost_yr"))
	{
		acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "ny");
		MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);
		for (ii=s0; ii<se; ii++)
		{
			ocp_nlp_cost_model_set(config, dims, in, ii, "y_ref", value);
		}
	}
	else if (!strcmp(field, "cost_yr_e"))
	{
		acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, N, "ny");
		MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);
		ocp_nlp_cost_model_set(config, dims, in, N, "y_ref", value);
	}
	else if (!strcmp(field, "init_x"))
	{
		if (nrhs!=6)
			MEX_SETTER_NO_STAGE_SUPPORT(fun_name, field)

		acados_size = (N+1) * nx;
		MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);
		for (ii=0; ii<=N; ii++)
		{
			ocp_nlp_out_set(config, dims, out, ii, "x", value+ii*nx);
		}
	}
	else if (!strcmp(field, "init_u"))
	{
		if (nrhs!=6)
			MEX_SETTER_NO_STAGE_SUPPORT(fun_name, field)

		acados_size = N*nu;
		MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);
		for (ii=0; ii<N; ii++)
		{
			ocp_nlp_out_set(config, dims, out, ii, "u", value+ii*nu);
		}
	}
	else if (!strcmp(field, "init_pi"))
	{
		if (nrhs!=6)
			MEX_SETTER_NO_STAGE_SUPPORT(fun_name, field)

		acados_size = N*nx;
		MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);
		for (ii=0; ii<N; ii++)
		{
			ocp_nlp_out_set(config, dims, out, ii, "pi", value+ii*nx);
		}
	}
	// TODO make setters for ALL numerical data
	else if (!strcmp(field, "p"))
	{
		external_function_param_casadi *ext_fun_param_ptr;

		// loop over number of external functions;
		int struct_size = mxGetNumberOfFields( C_ext_fun_pointers );
		for (ii=0; ii<struct_size; ii++)
		{
//			printf("\n%s\n", mxGetFieldNameByNumber( C_ext_fun_pointers, ii) );
			mex_field = mxGetFieldByNumber( C_ext_fun_pointers, 0, ii );
			ptr = (long long *) mxGetData( mex_field );
			int Nf_sum = 0;
			// loop over number of phases
			for (jj=0; jj<Nf; jj++)
			{
				ext_fun_param_ptr = (external_function_param_casadi *) ptr[jj];
				if (ext_fun_param_ptr!=0)
				{
					if (nrhs==6)
					{
						for (kk=0; kk<NN[jj]; kk++)
						{
							acados_size = (ext_fun_param_ptr+kk)->np;
							MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);
							(ext_fun_param_ptr+kk)->set_param(ext_fun_param_ptr+kk, value);
						}
					}
					else if (nrhs==7)
					{
						int stage = mxGetScalar( prhs[6] );
						if (stage>=Nf_sum & stage<Nf_sum+NN[jj])
						{
							acados_size = (ext_fun_param_ptr+stage)->np;
							MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);
							(ext_fun_param_ptr+stage)->set_param(ext_fun_param_ptr+stage, value);
						}
					}
				}
				Nf_sum += NN[jj];
			}
		}
	}
	else if (!strcmp(field, "nlp_solver_max_iter"))
	{
		acados_size = 1;
		MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);
		int nlp_solver_max_iter = (int) value[0];
		ocp_nlp_opts_set(config, opts, "max_iter", &nlp_solver_max_iter);
	}
	else
	{
		MEX_FIELD_NOT_SUPPORTED(fun_name, field);
	}

	return;
}

