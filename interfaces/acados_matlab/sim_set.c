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

	// solver
	ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "solver" ) );
	sim_solver *solver = (sim_solver *) ptr[0];
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
	else if(!strcmp(field, "xdot"))
		{
		double *xdot = mxGetPr( prhs[5] );
		sim_solver_set(solver, "xdot", xdot);
		}
	else if(!strcmp(field, "phi_guess"))
		{
		double *phi_guess = mxGetPr( prhs[5] );
		sim_solver_set(solver, "phi_guess", phi_guess);
		}
	else
		{
		mexPrintf("\nsim_set: field not supported: %s\n", field);
		return;
		}



	/* return */

	return;

	}


