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
#include "acados/utils/external_function_generic.h"
#include "acados_c/external_function_interface.h"
// mex
#include "mex.h"
#include "mex_macros.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    int ii, jj, kk;
    long long *ptr;
    mxArray *mex_field;
    char fun_name[20] = "ocp_destroy_ext_fun";
    char buffer [400]; // for error messages

    /* RHS */

    // model_struct
    char *cost_ext_fun_type;
    const mxArray *matlab_model = prhs[0];
    if (mxGetField( matlab_model, 0, "cost_ext_fun_type" )!=NULL)
        cost_ext_fun_type = mxArrayToString( mxGetField( matlab_model, 0, "cost_ext_fun_type" ) );

    char *cost_ext_fun_type_e;    
    if (mxGetField( matlab_model, 0, "cost_ext_fun_type_e" )!=NULL)
        cost_ext_fun_type_e = mxArrayToString( mxGetField( matlab_model, 0, "cost_ext_fun_type_e" ) );

    char *cost_ext_fun_type_0;    
    if (mxGetField( matlab_model, 0, "cost_ext_fun_type_0" )!=NULL)
        cost_ext_fun_type_0 = mxArrayToString( mxGetField( matlab_model, 0, "cost_ext_fun_type_0" ) );

    char *dyn_ext_fun_type;
    if (mxGetField( matlab_model, 0, "dyn_ext_fun_type" )!=NULL)
        dyn_ext_fun_type = mxArrayToString( mxGetField( matlab_model, 0, "dyn_ext_fun_type" ) );

    // dims
    ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dims" ) );
    ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];

    int N = dims->N;

    // XXX hard-code number and size of phases for now !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    int NN[] = {N, 1}; // number of phases, i.e. shooting nodes with same dimensions
    int Nf;

    //
    int struct_size = mxGetNumberOfFields( prhs[2] );
    for (ii=0; ii<struct_size; ii++)
    {
        //printf("\n%s\n", mxGetFieldNameByNumber( prhs[2], ii) );
        mex_field = mxGetFieldByNumber( prhs[2], 0, ii );
        ptr = (long long *) mxGetData( mex_field );
        Nf = mxGetN( mex_field );
        
        if (!strcmp(mxGetFieldNameByNumber(prhs[2], ii), "cost_ext_cost_fun") || 
            !strcmp(mxGetFieldNameByNumber(prhs[2], ii), "cost_ext_cost_fun_jac_hess"))
        {
            // TODO: what is Nf?? How should stage 0 be considered here???
            for (jj=0; jj<Nf; jj++)
            {
                    // external function param casadi
                    if ((jj == 0 && !strcmp(cost_ext_fun_type, "casadi")) || 
                        (jj == 1 && !strcmp(cost_ext_fun_type_e, "casadi")))
                    {
                        external_function_param_casadi *ext_fun_ptr = (external_function_param_casadi *) ptr[jj];
                        if (ext_fun_ptr!=0)
                        {
                            for (kk=0; kk<NN[jj]; kk++)
                            {
                                external_function_param_casadi_free(ext_fun_ptr+kk);
                            }
                            free(ext_fun_ptr);
                        }
                    }
                    // external function param generic
                    else if ((jj == 0 && !strcmp(cost_ext_fun_type, "generic")) || 
                             (jj == 1 && !strcmp(cost_ext_fun_type_e, "generic")))
                    {
                        external_function_param_generic *ext_fun_ptr = (external_function_param_generic *) ptr[jj];
                        if (ext_fun_ptr!=0)
                        {
                            for (kk=0; kk<NN[jj]; kk++)
                            {
                                external_function_param_generic_free(ext_fun_ptr+kk);
                            }
                            free(ext_fun_ptr);
                        }
                    }
                    else
                    {
                        MEX_FIELD_VALUE_NOT_SUPPORTED_SUGGEST(fun_name, "cost_ext_fun_type", cost_ext_fun_type, "casadi, generic");
                    }
            }
        }
        else if (!strcmp(mxGetFieldNameByNumber(prhs[2], ii), "dyn_disc_phi_fun") ||
                 !strcmp(mxGetFieldNameByNumber(prhs[2], ii), "dyn_disc_phi_fun_jac") ||
                 !strcmp(mxGetFieldNameByNumber(prhs[2], ii), "dyn_disc_phi_fun_jac_hess"))
        {
            for (jj=0; jj<Nf; jj++)
            {
                    // external function param casadi
                    if (!strcmp(dyn_ext_fun_type, "casadi"))
                    {
                        external_function_param_casadi *ext_fun_ptr = (external_function_param_casadi *) ptr[jj];
                        if (ext_fun_ptr!=0)
                        {
                            for (kk=0; kk<NN[jj]; kk++)
                            {
                                external_function_param_casadi_free(ext_fun_ptr+kk);
                            }
                            free(ext_fun_ptr);
                        }
                    }
                    // external function param generic
                    else if (!strcmp(dyn_ext_fun_type, "generic"))
                    {
                        external_function_param_generic *ext_fun_ptr = (external_function_param_generic *) ptr[jj];
                        if (ext_fun_ptr!=0)
                        {
                            for (kk=0; kk<NN[jj]; kk++)
                            {
                                external_function_param_generic_free(ext_fun_ptr+kk);
                            }
                            free(ext_fun_ptr);
                        }
                    }
                    else
                    {
                        MEX_FIELD_VALUE_NOT_SUPPORTED_SUGGEST(fun_name, "cost_ext_fun_type", cost_ext_fun_type, "casadi, generic");
                    }
            }
        }
        else
        {
            for (jj=0; jj<Nf; jj++)
            {
                external_function_param_casadi *ext_fun_ptr = (external_function_param_casadi *) ptr[jj];
                if (ext_fun_ptr!=0)
                {
                    for (kk=0; kk<NN[jj]; kk++)
                    {
                        external_function_param_casadi_free(ext_fun_ptr+kk);
                    }
                    free(ext_fun_ptr);
                }
            }
        }
    }

    return;

}


