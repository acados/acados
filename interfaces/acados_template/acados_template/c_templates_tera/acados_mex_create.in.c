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


// standard
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// acados
#include "acados/utils/print.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_solver_{{ model.name }}.h"

// mex
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    long long *l_ptr;
    int status = 0;

    // create solver
    status = acados_create();

    if (status)
    {
        mexPrintf("acados_create() returned status %d.\n", status);
    }
    // mexPrintf("acados_create -> success!\n");

    // get pointers to nlp solver related objects
    ocp_nlp_plan *nlp_plan = acados_get_nlp_plan();
    ocp_nlp_config *nlp_config = acados_get_nlp_config();
    ocp_nlp_dims *nlp_dims = acados_get_nlp_dims();
    ocp_nlp_in *nlp_in = acados_get_nlp_in();
    ocp_nlp_out *nlp_out = acados_get_nlp_out();
    ocp_nlp_solver *nlp_solver = acados_get_nlp_solver();
    void *nlp_opts = acados_get_nlp_opts();

    // mexPrintf("acados: got pointer to objectes!\n");

    // field names of output struct
    #define FIELDS_OCP 8
    #define FIELDS_EXT_FUN 7
    #define MAX_FIELDS 8
    char *fieldnames[MAX_FIELDS];

    for (int i = 0; i < MAX_FIELDS; i++)
    {
        fieldnames[i] = (char*) mxMalloc(50);
    }

    memcpy(fieldnames[0],"config",sizeof("config"));
    memcpy(fieldnames[1],"dims",sizeof("dims"));
    memcpy(fieldnames[2],"opts",sizeof("opts"));
    memcpy(fieldnames[3],"in",sizeof("in"));
    memcpy(fieldnames[4],"out",sizeof("out"));
    memcpy(fieldnames[5],"solver",sizeof("solver"));
    memcpy(fieldnames[6],"sens_out",sizeof("sens_out"));
    memcpy(fieldnames[7],"plan",sizeof("plan"));

    // create output struct - C_ocp
    plhs[0] = mxCreateStructMatrix(1, 1, 8, (const char **) fieldnames);

    // MEX: config, dims, opts, in, out, solver, sens_out, plan
    // plan
    mxArray *plan_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(plan_mat);
    l_ptr[0] = (long long) nlp_plan;
    mxSetField(plhs[0], 0, "plan", plan_mat);

    // config
    mxArray *config_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(config_mat);
    l_ptr[0] = (long long) nlp_config;
    mxSetField(plhs[0], 0, "config", config_mat);

    // dims
    mxArray *dims_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(dims_mat);
    l_ptr[0] = (long long) nlp_dims;
    mxSetField(plhs[0], 0, "dims", dims_mat);

    // opts
    mxArray *opts_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(opts_mat);
    l_ptr[0] = (long long) nlp_opts;
    mxSetField(plhs[0], 0, "opts", opts_mat);

    // in
    mxArray *in_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(in_mat);
    l_ptr[0] = (long long) nlp_in;
    mxSetField(plhs[0], 0, "in", in_mat);

    // out
    mxArray *out_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(out_mat);
    l_ptr[0] = (long long) nlp_out;
    mxSetField(plhs[0], 0, "out", out_mat);

    // solver
    mxArray *solver_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(solver_mat);
    l_ptr[0] = (long long) nlp_solver;
    mxSetField(plhs[0], 0, "solver", solver_mat);

    // TODO: sens_out not actually implemented in templates..
    // sens_out
    mxArray *sens_out_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(sens_out_mat);
    l_ptr[0] = (long long) 1;
    mxSetField(plhs[0], 0, "sens_out", sens_out_mat);

    /* store external function pointers */
    memcpy(fieldnames[0],"forw_vde",sizeof("forw_vde"));
    memcpy(fieldnames[1],"hess_vde",sizeof("hess_vde"));
    memcpy(fieldnames[2],"impl_dae_fun",sizeof("impl_dae_fun"));
    memcpy(fieldnames[3],"impl_dae_fun_jac_x_xdot_z",sizeof("impl_dae_fun_jac_x_xdot_z"));
    memcpy(fieldnames[4],"impl_dae_jac_x_xdot_u_z",sizeof("impl_dae_jac_x_xdot_u_z"));
    // constraints
    memcpy(fieldnames[5],"phi_constraint",sizeof("phi_constraint"));
    memcpy(fieldnames[6],"h_constraint",sizeof("h_constraint"));

    // create output struct - C_ocp_ext_fun
    plhs[1] = mxCreateStructMatrix(1, 1, FIELDS_EXT_FUN, (const char **) fieldnames);


    for (int i = 0; i < FIELDS_EXT_FUN; i++)
    {
        mxFree( fieldnames[i] );
    }

// dynamics
    mxArray *forw_vde_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    mxArray *hess_vde_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    mxArray *impl_dae_fun_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    mxArray *impl_dae_fun_jac_x_xdot_z_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    mxArray *impl_dae_jac_x_xdot_u_z_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);

{% if solver_options.integrator_type == "ERK" %}
    {# TODO: remove _casadi from these names.. #}
    l_ptr = mxGetData(forw_vde_mat);
    l_ptr[0] = (long long) forw_vde_casadi;
    // mexPrintf("\nforw vde %p\n", forw_vde_casadi);
{% if solver_options.hessian_approx == "EXACT" %}
    l_ptr = mxGetData(hess_vde_mat);
    l_ptr[0] = (long long) hess_vde_casadi;
{%- endif %}
{% elif solver_options.integrator_type == "IRK" %}
    // extern external_function_param_casadi * impl_dae_fun;
    l_ptr = mxGetData(impl_dae_fun_mat);
    l_ptr[0] = (long long) impl_dae_fun;
    // extern external_function_param_casadi * impl_dae_fun_jac_x_xdot_z;
    l_ptr = mxGetData(impl_dae_fun_jac_x_xdot_z_mat);
    l_ptr[0] = (long long) impl_dae_fun_jac_x_xdot_z;
    // extern external_function_param_casadi * impl_dae_jac_x_xdot_u_z;
    l_ptr = mxGetData(impl_dae_jac_x_xdot_u_z_mat);
    l_ptr[0] = (long long) impl_dae_jac_x_xdot_u_z;
{%- endif %}
    mxSetField(plhs[1], 0, "forw_vde", forw_vde_mat);
    mxSetField(plhs[1], 0, "hess_vde", hess_vde_mat);

    mxSetField(plhs[1], 0, "impl_dae_fun", impl_dae_fun_mat);
    mxSetField(plhs[1], 0, "impl_dae_fun_jac_x_xdot_z", impl_dae_fun_jac_x_xdot_z_mat);
    mxSetField(plhs[1], 0, "impl_dae_jac_x_xdot_u_z", impl_dae_jac_x_xdot_u_z_mat);

// constaints
    mxArray *phi_constraint_mat = mxCreateNumericMatrix(1, 2, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(phi_constraint_mat);
{%- if constraints.constr_type == "BGP" %}
    l_ptr[0] = (long long) phi_constraint;
{% endif %}
{% if constraints.constr_type_e == "BGP" %}
    l_ptr[1] = (long long) &phi_e_constraint;
{% endif %}
    mxSetField(plhs[1], 0, "phi_constraint", phi_constraint_mat);

    mxArray *h_constraint_mat  = mxCreateNumericMatrix(1, 2, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(h_constraint_mat);
{% if constraints.constr_type == "BGH" and dims.nh > 0 %}
    l_ptr[0] = (long long) h_constraint;
{% endif %}
{% if constraints.constr_type_e == "BGH" and dims.nh_e > 0 %}
    l_ptr[1] = (long long) &h_e_constraint;
{%- endif %}
    mxSetField(plhs[1], 0, "h_constraint", h_constraint_mat);


// {% if cost.cost_type == "NONLINEAR_LS" %}
//     mxArray *cost_y_fun_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
//     l_ptr = mxGetData(cost_y_fun_mat);
//     l_ptr[0] = (long long) cost_y_fun;
//     mxSetField(plhs[1], 0, "cost_y_fun", cost_y_fun_mat);
// {% endif %}
// {% if cost.cost_type_e == "NONLINEAR_LS" %}
//     mxArray *cost_y_e_fun_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
//     l_ptr = mxGetData(cost_y_e_fun_mat);
//     l_ptr[0] = (long long) cost_y_e_fun;
//     mxSetField(plhs[1], 0, "cost_y_e_fun", cost_y_e_fun_mat);
// {%- endif %}

    return;
}
