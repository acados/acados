/*
 * Copyright (c) The acados authors.
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
//#include "acados/utils/external_function_generic.h"
//#include "acados_c/external_function_interface.h"
// mex
#include "mex.h"



// casadi functions for the model
//#include "model.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    // field names of output struct
    char *fieldnames[24];
    fieldnames[0] = (char*)mxMalloc(50);
    fieldnames[1] = (char*)mxMalloc(50);
    fieldnames[2] = (char*)mxMalloc(50);
    fieldnames[3] = (char*)mxMalloc(50);
    fieldnames[4] = (char*)mxMalloc(50);
    fieldnames[5] = (char*)mxMalloc(50);
    fieldnames[6] = (char*)mxMalloc(50);
    fieldnames[7] = (char*)mxMalloc(50);
    fieldnames[8] = (char*)mxMalloc(50);
    fieldnames[9] = (char*)mxMalloc(50);
    fieldnames[10] = (char*)mxMalloc(50);
    fieldnames[11] = (char*)mxMalloc(50);
    fieldnames[12] = (char*)mxMalloc(50);
    fieldnames[13] = (char*)mxMalloc(50);
    fieldnames[14] = (char*)mxMalloc(50);
    fieldnames[15] = (char*)mxMalloc(50);
    fieldnames[16] = (char*)mxMalloc(50);
    fieldnames[17] = (char*)mxMalloc(50);
    fieldnames[18] = (char*)mxMalloc(50);
    fieldnames[19] = (char*)mxMalloc(50);
    fieldnames[20] = (char*)mxMalloc(50);
    fieldnames[21] = (char*)mxMalloc(50);
    fieldnames[22] = (char*)mxMalloc(50);
    fieldnames[23] = (char*)mxMalloc(50);

    memcpy(fieldnames[0],"dyn_expl_ode_fun",sizeof("dyn_expl_ode_fun"));
    memcpy(fieldnames[1],"dyn_expl_vde_forw",sizeof("dyn_expl_vde_forw"));
    memcpy(fieldnames[2],"dyn_expl_vde_adj",sizeof("dyn_expl_vde_adj"));
    memcpy(fieldnames[3],"dyn_expl_ode_hess",sizeof("dyn_expl_ode_hess"));
    memcpy(fieldnames[4],"dyn_impl_dae_fun",sizeof("dyn_impl_dae_fun"));
    memcpy(fieldnames[5],"dyn_impl_dae_fun_jac_x_xdot_z",sizeof("dyn_impl_dae_fun_jac_x_xdot_z"));
    memcpy(fieldnames[6],"dyn_impl_dae_jac_x_xdot_u_z",sizeof("dyn_impl_dae_jac_x_xdot_u_z"));
    memcpy(fieldnames[7],"dyn_impl_dae_hess",sizeof("dyn_impl_dae_hess"));
    memcpy(fieldnames[8],"dyn_gnsf_f_lo_fun_jac_x1k1uz",sizeof("dyn_gnsf_f_lo_fun_jac_x1k1uz"));
    memcpy(fieldnames[9],"dyn_gnsf_get_matrices_fun",sizeof("dyn_gnsf_get_matrices_fun"));
    memcpy(fieldnames[10],"dyn_gnsf_phi_fun",sizeof("dyn_gnsf_phi_fun"));
    memcpy(fieldnames[11],"dyn_gnsf_phi_fun_jac_y",sizeof("dyn_gnsf_phi_fun_jac_y"));
    memcpy(fieldnames[12],"dyn_gnsf_phi_jac_y_uhat",sizeof("dyn_gnsf_phi_jac_y_uhat"));
    memcpy(fieldnames[13],"dyn_disc_phi_fun",sizeof("dyn_disc_phi_fun"));
    memcpy(fieldnames[14],"dyn_disc_phi_fun_jac",sizeof("dyn_disc_phi_fun_jac"));
    memcpy(fieldnames[15],"dyn_disc_phi_fun_jac_hess",sizeof("dyn_disc_phi_fun_jac_hess"));
    memcpy(fieldnames[16],"constr_h_fun",sizeof("constr_h_fun"));
    memcpy(fieldnames[17],"constr_h_fun_jac_uxt_zt",sizeof("constr_h_fun_jac_uxt_zt"));
    memcpy(fieldnames[18],"constr_h_fun_jac_uxt_zt_hess",sizeof("constr_h_fun_jac_uxt_zt_hess"));
    memcpy(fieldnames[19],"cost_y_fun",sizeof("cost_y_fun"));
    memcpy(fieldnames[20],"cost_y_fun_jac_ut_xt",sizeof("cost_y_fun_jac_ut_xt"));
    memcpy(fieldnames[21],"cost_y_hess",sizeof("cost_y_hess"));
    memcpy(fieldnames[22],"cost_ext_cost_fun",sizeof("cost_ext_cost_fun"));
    memcpy(fieldnames[23],"cost_ext_cost_fun_jac_hess",sizeof("cost_ext_cost_fun_jac_hess"));

    // create output struct
    plhs[0] = mxCreateStructMatrix(1, 1, 24, (const char **) fieldnames);

    mxFree( fieldnames[0] );
    mxFree( fieldnames[1] );
    mxFree( fieldnames[2] );
    mxFree( fieldnames[3] );
    mxFree( fieldnames[4] );
    mxFree( fieldnames[5] );
    mxFree( fieldnames[6] );
    mxFree( fieldnames[7] );
    mxFree( fieldnames[8] );
    mxFree( fieldnames[9] );
    mxFree( fieldnames[10] );
    mxFree( fieldnames[11] );
    mxFree( fieldnames[12] );
    mxFree( fieldnames[13] );
    mxFree( fieldnames[14] );
    mxFree( fieldnames[15] );
    mxFree( fieldnames[16] );
    mxFree( fieldnames[17] );
    mxFree( fieldnames[18] );
    mxFree( fieldnames[19] );
    mxFree( fieldnames[20] );
    mxFree( fieldnames[21] );
    mxFree( fieldnames[22] );
    mxFree( fieldnames[23] );

    // populate struct with empty vectors with number of phases length
    int Nf = 1;

    mxSetField(plhs[0], 0, "dyn_expl_ode_fun", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "dyn_expl_vde_forw", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "dyn_expl_vde_adj", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "dyn_expl_ode_hess", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "dyn_impl_dae_fun", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "dyn_impl_dae_fun_jac_x_xdot_z", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "dyn_impl_dae_jac_x_xdot_u_z", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "dyn_impl_dae_hess", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "dyn_gnsf_f_lo_fun_jac_x1k1uz", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "dyn_gnsf_get_matrices_fun", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "dyn_gnsf_phi_fun", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "dyn_gnsf_phi_fun_jac_y", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "dyn_gnsf_phi_jac_y_uhat", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "dyn_disc_phi_fun", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "dyn_disc_phi_fun_jac", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "dyn_disc_phi_fun_jac_hess", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "constr_h_fun", mxCreateNumericMatrix(1, Nf+1, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "constr_h_fun_jac_uxt_zt", mxCreateNumericMatrix(1, Nf+1, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "constr_h_fun_jac_uxt_zt_hess", mxCreateNumericMatrix(1, Nf+1, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "cost_y_fun", mxCreateNumericMatrix(1, Nf+1, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "cost_y_fun_jac_ut_xt", mxCreateNumericMatrix(1, Nf+1, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "cost_y_hess", mxCreateNumericMatrix(1, Nf+1, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "cost_ext_cost_fun", mxCreateNumericMatrix(1, Nf+1, mxINT64_CLASS, mxREAL));
    mxSetField(plhs[0], 0, "cost_ext_cost_fun_jac_hess", mxCreateNumericMatrix(1, Nf+1, mxINT64_CLASS, mxREAL));

    return;

}



