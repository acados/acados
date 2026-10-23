/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */

#define S_FUNCTION_NAME   acados_sim_solver_sfunction_{{ name }}
#define S_FUNCTION_LEVEL  2

#define MDL_START

// acados
// #include "acados/utils/print.h"
#include "acados_c/sim_interface.h"
#include "acados_c/external_function_interface.h"

// example specific
#include "{{ model.name }}_model/{{ model.name }}_model.h"
#include "acados_sim_solver_{{ name }}.h"

#include "simstruc.h"

#define SAMPLINGTIME {{ solver_options.Tsim }}


typedef struct {
    {{ name }}_sim_solver_capsule *capsule;
    double* buffer;
} AcadosSimData;


static void mdlInitializeSizes (SimStruct *S)
{
    // specify the number of continuous and discrete states
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    {# compute number of input ports #}
    {%- set n_inputs = 0 %}
    {%- if simulink_opts.inputs.x %}{%- set n_inputs = n_inputs + 1 %}{%- endif %}
    {%- if dims.nu > 0 and simulink_opts.inputs.u %}{%- set n_inputs = n_inputs + 1 %}{%- endif %}
    {%- if dims.np > 0 and simulink_opts.inputs.p %}{%- set n_inputs = n_inputs + 1 %}{%- endif %}
    {%- if simulink_opts.inputs.xdot %}{%- set n_inputs = n_inputs + 1 %}{%- endif %}
    {%- if dims.nz > 0 and simulink_opts.inputs.z %}{%- set n_inputs = n_inputs + 1 %}{%- endif %}
    {%- if simulink_opts.inputs.seed_adj %}{%- set n_inputs = n_inputs + 1 %}{%- endif %}
    {%- if simulink_opts.inputs.T %}{%- set n_inputs = n_inputs + 1 %}{%- endif %}
    {%- if simulink_opts.inputs.t0 %}{%- set n_inputs = n_inputs + 1 %}{%- endif %}

    {%- set n_outputs = 0 %}
    {%- if simulink_opts.outputs.x %}{%- set n_outputs = n_outputs + 1 %}{%- endif %}
    {%- if dims.nz > 0 and simulink_opts.outputs.z %}{%- set n_outputs = n_outputs + 1 %}{%- endif %}
    {%- if simulink_opts.outputs.S_forw %}{%- set n_outputs = n_outputs + 1 %}{%- endif %}
    {%- if simulink_opts.outputs.S_adj %}{%- set n_outputs = n_outputs + 1 %}{%- endif %}
    {%- if simulink_opts.outputs.time_tot %}{%- set n_outputs = n_outputs + 1 %}{%- endif %}

    // specify the number of input ports
    if ( !ssSetNumInputPorts(S, {{ n_inputs }}) )
        return;

    // specify the number of output ports
    if ( !ssSetNumOutputPorts(S, {{ n_outputs }}) )
        return;

    // specify dimension information for the input ports
    {%- set i_input = -1 %}
    {%- if simulink_opts.inputs.x %}{%- set i_input = i_input + 1 %}
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nx }});
    {%- endif %}
    {%- if dims.nu > 0 and simulink_opts.inputs.u %}{%- set i_input = i_input + 1 %}
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nu }});
    {%- endif %}
    {%- if dims.np > 0 and simulink_opts.inputs.p %}{%- set i_input = i_input + 1 %}
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.np }});
    {%- endif %}
    {%- if simulink_opts.inputs.xdot %}{%- set i_input = i_input + 1 %}
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nx }});
    {%- endif %}
    {%- if dims.nz > 0 and simulink_opts.inputs.z %}{%- set i_input = i_input + 1 %}
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nz }});
    {%- endif %}
    {%- if simulink_opts.inputs.seed_adj %}{%- set i_input = i_input + 1 %}
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nx }});
    {%- endif %}
    {%- if simulink_opts.inputs.T %}{%- set i_input = i_input + 1 %}
    ssSetInputPortVectorDimension(S, {{ i_input }}, 1);
    {%- endif %}
    {%- if simulink_opts.inputs.t0 %}{%- set i_input = i_input + 1 %}
    ssSetInputPortVectorDimension(S, {{ i_input }}, 1);
    {%- endif %}

    {%- set i_output = -1 %}
    {%- if simulink_opts.outputs.x %}{%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, {{ dims.nx }});
    {%- endif %}
    {%- if dims.nz > 0 and simulink_opts.outputs.z %}{%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, {{ dims.nz }});
    {%- endif %}
    {%- if simulink_opts.outputs.S_forw %}{%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, {{ dims.nx * (dims.nx + dims.nu) }});
    {%- endif %}
    {%- if simulink_opts.outputs.S_adj %}{%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, {{ dims.nx + dims.nu }});
    {%- endif %}
    {%- if simulink_opts.outputs.time_tot %}{%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, 1);
    {%- endif %}

    // specify the direct feedthrough status
    // should be set to 1 for all inputs used in mdlOutputs
    {%- for i in range(end=n_inputs) %}
    ssSetInputPortDirectFeedThrough(S, {{ i }}, 1);
    {%- endfor %}

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
    AcadosSimData *sim_data = malloc(sizeof(*sim_data));

    // capsule
    sim_data->capsule = {{ name }}_acados_sim_solver_create_capsule();
    {{ name }}_acados_sim_create(sim_data->capsule);

    // local buffer
    {% set input_sizes = [dims.nx, dims.nu, dims.np, dims.nz, dims.nx * (dims.nx + dims.nu), 1] %}
    {%- set buffer_size =  input_sizes | sort | last %}
    sim_data->buffer = malloc({{ buffer_size}} * sizeof(double));

    ssSetUserData(S, (void*) sim_data);
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    AcadosSimData *sim_data = ssGetUserData(S);
    {{ name }}_sim_solver_capsule *capsule = sim_data->capsule;
    double* buffer = sim_data->buffer;

    sim_config *acados_sim_config = {{ name }}_acados_get_sim_config(capsule);
    sim_in *acados_sim_in = {{ name }}_acados_get_sim_in(capsule);
    sim_out *acados_sim_out = {{ name }}_acados_get_sim_out(capsule);
    void *acados_sim_dims = {{ name }}_acados_get_sim_dims(capsule);
    sim_solver *acados_sim_solver = {{ name }}_acados_get_sim_solver(capsule);

    InputRealPtrsType in_sign;

    /* go through inputs */
    {%- set i_input = -1 %}
{% if simulink_opts.inputs.x %}
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int i = 0; i < {{ dims.nx }}; i++) buffer[i] = (double)(*in_sign[i]);
    sim_in_set(acados_sim_config, acados_sim_dims, acados_sim_in, "x", buffer);
{% endif %}
{% if dims.nu > 0 and simulink_opts.inputs.u %}
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int i = 0; i < {{ dims.nu }}; i++) buffer[i] = (double)(*in_sign[i]);
    sim_in_set(acados_sim_config, acados_sim_dims, acados_sim_in, "u", buffer);
{% endif %}
{% if dims.np > 0 and simulink_opts.inputs.p %}
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int i = 0; i < {{ dims.np }}; i++) buffer[i] = (double)(*in_sign[i]);
    {{ name }}_acados_sim_update_params(capsule, buffer, {{ dims.np }});
{% endif %}
{% if simulink_opts.inputs.xdot %}
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int i = 0; i < {{ dims.nx }}; i++) buffer[i] = (double)(*in_sign[i]);
    sim_solver_set(acados_sim_solver, "xdot", buffer);
{% endif %}
{% if dims.nz > 0 and simulink_opts.inputs.z %}
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int i = 0; i < {{ dims.nz }}; i++) buffer[i] = (double)(*in_sign[i]);
    sim_solver_set(acados_sim_solver, "z", buffer);
{% endif %}
{% if simulink_opts.inputs.seed_adj %}
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int i = 0; i < {{ dims.nx }}; i++) buffer[i] = (double)(*in_sign[i]);
    sim_in_set(acados_sim_config, acados_sim_dims, acados_sim_in, "seed_adj", buffer);
{% endif %}
{% if simulink_opts.inputs.T %}
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    buffer[0] = (double)(*in_sign[0]);
    sim_in_set(acados_sim_config, acados_sim_dims, acados_sim_in, "T", buffer);
{% endif %}
{% if simulink_opts.inputs.t0 %}
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    buffer[0] = (double)(*in_sign[0]);
    sim_in_set(acados_sim_config, acados_sim_dims, acados_sim_in, "t0", buffer);
{% endif %}


    /* call solver */
    int acados_status = {{ name }}_acados_sim_solve(capsule);


    /* set outputs */
    {%- set i_output = -1 %}
{% if simulink_opts.outputs.x %}
    {%- set i_output = i_output + 1 %}
    double *out_x = ssGetOutputPortRealSignal(S, {{ i_output }});
    sim_out_get(acados_sim_config, acados_sim_dims, acados_sim_out, "x", (void *) out_x);
{% endif %}
{% if dims.nz > 0 and simulink_opts.outputs.z %}
    {%- set i_output = i_output + 1 %}
    double *out_z = ssGetOutputPortRealSignal(S, {{ i_output }});
    sim_out_get(acados_sim_config, acados_sim_dims, acados_sim_out, "z", (void *) out_z);
{% endif %}
{% if simulink_opts.outputs.S_forw %}
    {%- set i_output = i_output + 1 %}
    double *out_S_forw = ssGetOutputPortRealSignal(S, {{ i_output }});
    sim_out_get(acados_sim_config, acados_sim_dims, acados_sim_out, "S_forw", (void *) out_S_forw);
{% endif %}
{% if simulink_opts.outputs.S_adj %}
    {%- set i_output = i_output + 1 %}
    double *out_S_adj = ssGetOutputPortRealSignal(S, {{ i_output }});
    sim_out_get(acados_sim_config, acados_sim_dims, acados_sim_out, "S_adj", (void *) out_S_adj);
{% endif %}
{% if simulink_opts.outputs.time_tot %}
    {%- set i_output = i_output + 1 %}
    double *out_time_tot = ssGetOutputPortRealSignal(S, {{ i_output }});
    sim_out_get(acados_sim_config, acados_sim_dims, acados_sim_out, "time_tot", (void *) out_time_tot);
{% endif %}

}


static void mdlTerminate(SimStruct *S)
{
    AcadosSimData *sim_data = ssGetUserData(S);
    {{ name }}_sim_solver_capsule *capsule = sim_data->capsule;
    double* buffer = sim_data->buffer;

    {{ name }}_acados_sim_free(capsule);
    {{ name }}_acados_sim_solver_free_capsule(capsule);
    free(buffer);
    free(sim_data);
}


#ifdef  MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
