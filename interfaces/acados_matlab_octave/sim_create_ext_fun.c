/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */

// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//    mexPrintf("\nin sim_expl_ext_fun_create\n");

    /* RHS */

    /* LHS */
    // field names of output struct
    char *fieldnames[15];
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
    fieldnames[13] = (char*)mxMalloc(50); // dyn_expl_vde_forw_p
    fieldnames[14] = (char*)mxMalloc(50); // dyn_impl_dae_jac_p

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
    memcpy(fieldnames[13],"dyn_expl_vde_forw_p",sizeof("dyn_expl_vde_forw_p")); // param-direction VDE
    memcpy(fieldnames[14],"dyn_impl_dae_jac_p",sizeof("dyn_impl_dae_jac_p"));

    // create output struct
    plhs[0] = mxCreateStructMatrix(1, 1, 15, (const char **) fieldnames);

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
    mxFree( fieldnames[13] ); // dyn_expl_vde_forw_p
    mxFree( fieldnames[14] );

//    mxSetField(plhs[0], 0, "dyn_impl_dae_fun", NULL );

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
    mxSetField(plhs[0], 0, "dyn_expl_vde_forw_p", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL)); // param-direction VDE
    mxSetField(plhs[0], 0, "dyn_impl_dae_jac_p", mxCreateNumericMatrix(1, Nf, mxINT64_CLASS, mxREAL));

    return;

}

