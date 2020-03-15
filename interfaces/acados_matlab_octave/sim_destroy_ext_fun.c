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
#include "acados/utils/external_function_generic.h"
#include "acados_c/external_function_interface.h"
// mex
#include "mex.h"
#include "mex_macros.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//    mexPrintf("\nin sim_ext_fun_destroy\n");

    int ii;
    long long *ptr;
//    mxArray *mex_field;
    char fun_name[20] = "sim_destroy_ext_fun";
    char buffer [400]; // for error messages

    /* RHS */

	// model_struct
	char *ext_fun_type;
    const mxArray *matlab_model = prhs[0];
    if (mxGetField( matlab_model, 0, "ext_fun_type" )!=NULL)
        ext_fun_type = mxArrayToString( mxGetField( matlab_model, 0, "ext_fun_type" ) );

    // C_sim_ext_fun

    int struct_size = mxGetNumberOfFields( prhs[1] );
    for (ii=0; ii<struct_size; ii++)
    {
//        printf("\n%s\n", mxGetFieldNameByNumber( prhs[1], ii) );
        ptr = (long long *) mxGetData( mxGetFieldByNumber( prhs[1], 0, ii ) );

		// external function param casadi
		if (!strcmp(ext_fun_type, "casadi"))
		{
			external_function_param_casadi *ext_fun_ptr = (external_function_param_casadi *) ptr[0];
			if (ext_fun_ptr!=0)
			{
				external_function_param_casadi_free(ext_fun_ptr);
				free(ext_fun_ptr);
			}
		}
		// external function param generic
		else if (!strcmp(ext_fun_type, "generic"))
		{
			external_function_param_generic *ext_fun_ptr = (external_function_param_generic *) ptr[0];
			if (ext_fun_ptr!=0)
			{
				external_function_param_generic_free(ext_fun_ptr);
				free(ext_fun_ptr);
			}
		}
		else
		{
			MEX_FIELD_VALUE_NOT_SUPPORTED_SUGGEST(fun_name, "ext_fun_type", ext_fun_type, "casadi, generic");
		}
    }

    return;
}

