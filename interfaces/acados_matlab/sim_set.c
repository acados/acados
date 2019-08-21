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
#include "acados/sim/sim_common.h"
#include "acados_c/sim_interface.h"
#include "acados/utils/external_function_generic.h"
#include "acados_c/external_function_interface.h"

// mex
#include "mex.h"
#include "mex_macros.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

	int ii;
	int acados_size, tmp;
	long long *ptr;
	char fun_name[50] = "sim_set";
	char buffer [100];

	/* RHS */
	// char *method = mxArrayToString( mxGetField( prhs[1], 0, "method" ) );

	// C_sim
	const mxArray *C_ext_fun_pointers = prhs[3];

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
	double *value = mxGetPr( prhs[5] );

	int matlab_size = (int) mxGetNumberOfElements( prhs[5] );

	// check dimension, set value
	if (!strcmp(field, "T"))
	{
		acados_size = 1;
		MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);
		sim_in_set(config, dims, in, field, value);
	}
	else if (!strcmp(field, "x"))
	{
		sim_dims_get(config, dims, "nx", &acados_size);
		MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);
		sim_in_set(config, dims, in, field, value);
	}
	else if (!strcmp(field, "u"))
	{
		sim_dims_get(config, dims, "nu", &acados_size);
		MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);
		sim_in_set(config, dims, in, field, value);
	}
	else if (!strcmp(field, "p"))
	{
		double *p = mxGetPr( prhs[5] );
		external_function_param_casadi *ext_fun_param_ptr;

		acados_size = ext_fun_param_ptr->np;
		MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);

		// loop over number of external functions;
		int struct_size = mxGetNumberOfFields( C_ext_fun_pointers );
		for(ii=0; ii<struct_size; ii++)
		{
//			printf("\n%s\n", mxGetFieldNameByNumber( C_ext_fun_pointers, ii) );
			ptr = (long long *) mxGetData( mxGetFieldByNumber( C_ext_fun_pointers, 0, ii ) );
			if (ptr[0]!=0)
			{
				ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
				ext_fun_param_ptr->set_param(ext_fun_param_ptr, p);
			}
		}
	}
	else if (!strcmp(field, "seed_adj"))
	{
		double *seed_adj = value;
		sim_dims_get(config, dims, "nx", &acados_size);
		// TODO(oj): in C, the backward seed is of dimension nx+nu, I think it should only be nx.
		// sim_dims_get(config, dims, "nu", &tmp);
		// acados_size += tmp;
		MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);
		sim_in_set(config, dims, in, "seed_adj", seed_adj);
	}
	else
	{
		sprintf(buffer, "sim_set: field not supported: %s\n", field);
		mexErrMsgTxt(buffer);
	}

	return;
}

