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

#define S_FUNCTION_NAME   acados_solver_sfunction_{{model.name}}
#define S_FUNCTION_LEVEL  2

#define MDL_START

// acados
#include "acados/utils/print.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

// TODO(oj) remove, when setters for Cyt,idxb available
#include "acados/ocp_nlp/ocp_nlp_constraints_bgh.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// example specific
#include "{{ model.name }}_model/{{ model.name }}_model.h"
#include "acados_solver_{{ model.name }}.h"

#include "simstruc.h"

#define SAMPLINGTIME -1
// ** global data **
ocp_nlp_in * nlp_in;
ocp_nlp_out * nlp_out;
ocp_nlp_solver * nlp_solver;
void * nlp_opts;
ocp_nlp_plan * nlp_solver_plan;
ocp_nlp_config * nlp_config;
ocp_nlp_dims * nlp_dims;
{% if solver_config.integrator_type == 'ERK' %}
{% if dims.np < 1 %}
external_function_casadi * forw_vde_casadi;
{% else %}
external_function_param_casadi * forw_vde_casadi;
{% endif %}
{% if solver_config.hessian_approx == 'EXACT' %} 
{% if dims.np < 1 %}
external_function_casadi * hess_vde_casadi;
{% else %}
external_function_param_casadi * hess_vde_casadi;
{% endif %}
{% endif %}
{% elif solver_config.integrator_type == 'IRK' %}
{% if dims.np < 1 %}
external_function_casadi * impl_dae_fun;
external_function_casadi * impl_dae_fun_jac_x_xdot_z;
external_function_casadi * impl_dae_jac_x_xdot_u_z;
{% else %}
external_function_param_casadi * impl_dae_fun;
external_function_param_casadi * impl_dae_fun_jac_x_xdot_z;
external_function_param_casadi * impl_dae_jac_x_xdot_u_z;
{% endif %}
{% endif %}
{% if dims.npd > 0 %}
external_function_casadi * p_constraint;
{% endif %}
{% if dims.npd_e > 0 %}
external_function_casadi p_e_constraint;
{% endif %}
{% if dims.nh > 0 %}
external_function_casadi * h_constraint;
{% endif %}
{% if dims.nh_e > 0 %}
external_function_casadi h_e_constraint;
{% endif %}


static void mdlInitializeSizes (SimStruct *S)
{
    // specify the number of continuous and discrete states 
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    // specify the number of input ports 
    {% if dims.np > 0 %}
    if ( !ssSetNumInputPorts(S, 4) )
    {% else %}
    if ( !ssSetNumInputPorts(S, 3) )
    {% endif %}
        return;

    // specify the number of output ports 
    if ( !ssSetNumOutputPorts(S, 5) )
        return;

    // specify dimension information for the input ports 
    ssSetInputPortVectorDimension(S, 0, {{ dims.nx }});
    ssSetInputPortVectorDimension(S, 1, {{ dims.ny }});
    ssSetInputPortVectorDimension(S, 2, {{ dims.ny_e }});
    {% if dims.np > 0 %}
    ssSetInputPortVectorDimension(S, 3, {{ dims.np }});
    {% endif %}

    // specify dimension information for the output ports 
    ssSetOutputPortVectorDimension(S, 0, {{ dims.nu }} ); // optimal input
    ssSetOutputPortVectorDimension(S, 1, 1 );                // solver status
    ssSetOutputPortVectorDimension(S, 2, 1 );                // KKT residuals
    ssSetOutputPortVectorDimension(S, 3, {{ dims.nx }} ); // first state
    ssSetOutputPortVectorDimension(S, 4, 1); // computation times

    // specify the direct feedthrough status 
    ssSetInputPortDirectFeedThrough(S, 0, 1); // current state x0
    ssSetInputPortDirectFeedThrough(S, 1, 1); // y_ref
    ssSetInputPortDirectFeedThrough(S, 2, 1); // y_ref_e
    {% if dims.np > 0 %}
    ssSetInputPortDirectFeedThrough(S, 3, 1); // parameter
    {% endif %}

    // one sample time 
    ssSetNumSampleTimes(S, 1);
    }


#if defined(MATLAB_MEX_FILE)

#define MDL_SET_INPUT_PORT_DIMENSION_INFO
#define MDL_SET_OUTPUT_PORT_DIMENSION_INFO

static void mdlSetInputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    if ( !ssSetInputPortDimensionInfo(S, port, dimsInfo) )
         return;
}

static void mdlSetOutputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    if ( !ssSetOutputPortDimensionInfo(S, port, dimsInfo) )
         return;
}

    #endif /* MATLAB_MEX_FILE */


static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, SAMPLINGTIME);
    ssSetOffsetTime(S, 0, 0.0);
}


static void mdlStart(SimStruct *S)
{
    acados_create();
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    // get input signals
    InputRealPtrsType in_x0_sign;
    InputRealPtrsType in_y_ref_sign;
    InputRealPtrsType in_y_ref_e_sign;
    {% if dims.np > 0 %}
    InputRealPtrsType in_p_sign;
    {% endif %}
    
    // local buffers
    real_t in_x0[{{ dims.nx }}];
    real_t in_y_ref[{{ dims.ny }}];
    real_t in_y_ref_e[{{ dims.ny_e }}];
    {% if dims.np > 0 %}
    real_t in_p[{{ dims.np }}];
    {% endif %}

    in_x0_sign = ssGetInputPortRealSignalPtrs(S, 0);
    in_y_ref_sign = ssGetInputPortRealSignalPtrs(S, 1);
    in_y_ref_e_sign = ssGetInputPortRealSignalPtrs(S, 2);
    {% if dims.np > 0 %}
    in_p_sign = ssGetInputPortRealSignalPtrs(S, 3);
    {% endif %}

    // copy signals into local buffers
    for (int i = 0; i < {{ dims.nx }}; i++) in_x0[i] = (double)(*in_x0_sign[i]);
    for (int i = 0; i < {{ dims.ny }}; i++) in_y_ref[i] = (double)(*in_y_ref_sign[i]);
    for (int i = 0; i < {{ dims.ny_e }}; i++) in_y_ref_e[i] = (double)(*in_y_ref_e_sign[i]);
    {% if dims.np > 0 %}
    for (int i = 0; i < {{ dims.np }}; i++) in_p[i] = (double)(*in_p_sign[i]);
    {% endif %}

    // for (int i = 0; i < 4; i++) ssPrintf("x0[%d] = %f\n", i, in_x0[i]);
    // ssPrintf("\n");

    // set initial condition
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", in_x0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", in_x0);

    // update reference
    for (int ii = 0; ii < {{dims.N}}; ii++)
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, ii, "yref", (void *) in_y_ref);

    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, {{dims.N}}, "yref", (void *) in_y_ref_e);

    // update value of parameters
    {% if dims.np > 0%}
    {% if solver_config.integrator_type == 'IRK' %}
    for (int ii = 0; ii < {{dims.N}}; ii++) {
    impl_dae_fun[ii].set_param(impl_dae_fun+ii, in_p);
    impl_dae_fun_jac_x_xdot_z[ii].set_param(impl_dae_fun_jac_x_xdot_z+ii, in_p);
    impl_dae_jac_x_xdot_u_z[ii].set_param(impl_dae_jac_x_xdot_u_z+ii, in_p);
    }
    {% else %}
    for (int ii = 0; ii < {{dims.N}}; ii++) {
    forw_vde_casadi[ii].set_param(forw_vde_casadi+ii, in_p);
    }
    {% endif %}
    {% endif %}
    
    // assign pointers to output signals 
    real_t *out_u0, *out_status, *out_KKT_res, *out_x1, *out_cpu_time;

    out_u0          = ssGetOutputPortRealSignal(S, 0);
    out_status      = ssGetOutputPortRealSignal(S, 1);
    out_KKT_res     = ssGetOutputPortRealSignal(S, 2);
    out_x1          = ssGetOutputPortRealSignal(S, 3);
    out_cpu_time    = ssGetOutputPortRealSignal(S, 4);
    
    // call acados_solve()
    int acados_status = acados_solve();

    *out_status = (real_t) acados_status;
    *out_KKT_res = (real_t) nlp_out->inf_norm_res;
    *out_cpu_time = (real_t) nlp_out->total_time;
    
    // get solution
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 0, "u", (void *) out_u0);

    // get next state
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 1, "x", (void *) out_x1);

}

static void mdlTerminate(SimStruct *S)
{
    acados_free();
}


#ifdef  MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
