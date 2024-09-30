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

#define S_FUNCTION_NAME acados_solver_sfunction_{{ name }}
#define S_FUNCTION_LEVEL 2

#define MDL_START

// acados
// #include "acados/utils/print.h"
#include "acados_c/sim_interface.h"
#include "acados_c/external_function_interface.h"

// example specific
#include "{{ name }}_model/{{ name }}_model.h"
#include "acados_solver_{{ name }}.h"


{%- if not solver_options.custom_update_filename %}
    {%- set custom_update_filename = "" %}
{% else %}
    {%- set custom_update_filename = solver_options.custom_update_filename %}
{%- endif %}

{%- if not solver_options.custom_update_header_filename %}
    {%- set custom_update_header_filename = "" %}
{% else %}
    {%- set custom_update_header_filename = solver_options.custom_update_header_filename %}
{%- endif %}
{%- if custom_update_header_filename != "" %}
#include "{{ custom_update_header_filename }}"
{%- endif %}

#include "simstruc.h"

{% if simulink_opts.samplingtime == "t0" -%}
#define SAMPLINGTIME {{ solver_options.time_steps[0] }}
{%- elif simulink_opts.samplingtime == "-1" -%}
#define SAMPLINGTIME -1
{%- else -%}
  {{ throw(message = "simulink_opts.samplingtime must be '-1' or 't0', got val") }}
{%- endif %}

{%- set ns_total = dims.ns_0 + dims.ns_e + (solver_options.N_horizon - 1) * dims.ns %}

static void mdlInitializeSizes (SimStruct *S)
{
    // specify the number of continuous and discrete states
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    int N = {{ model.name | upper }}_N;

  {%- for key, val in simulink_opts.inputs -%}
    {%- if val != 0 and val != 1 -%}
      {{ throw(message = "simulink_opts.inputs must be 0 or 1, got val") }}
    {%- endif -%}
  {%- endfor -%}

  {#- compute number of input ports #}
  {%- set n_inputs = 0 -%}
  {%- if dims.nbx_0 > 0 and simulink_opts.inputs.lbx_0 -%}  {#- lbx_0 #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.nbx_0 > 0 and simulink_opts.inputs.ubx_0 -%}  {#- ubx_0 #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.np > 0 and simulink_opts.inputs.parameter_traj -%}  {#- parameter_traj #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.ny_0 > 0 and simulink_opts.inputs.y_ref_0 -%}  {#- y_ref_0 -#}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.ny > 0 and solver_options.N_horizon > 1 and simulink_opts.inputs.y_ref -%}  {#- y_ref -#}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.ny_e > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.y_ref_e -%}  {#- y_ref_e #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.nbx > 0 and solver_options.N_horizon > 1 and simulink_opts.inputs.lbx -%}  {#- lbx #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.nbx > 0 and solver_options.N_horizon > 1 and simulink_opts.inputs.ubx -%}  {#- ubx #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.nbx_e > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.lbx_e -%}  {#- lbx_e #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.nbx_e > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.ubx_e -%}  {#- ubx_e #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.nbu > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.lbu -%}  {#- lbu #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.nbu > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.ubu -%}  {#- ubu #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.ng > 0 and simulink_opts.inputs.lg -%}  {#- lg #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.ng > 0 and simulink_opts.inputs.ug -%}  {#- ug #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.nh > 0 and simulink_opts.inputs.lh -%}  {#- lh #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.nh > 0 and simulink_opts.inputs.uh -%}  {#- uh #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.nh_0 > 0 and simulink_opts.inputs.lh_0 -%}  {#- lh_0 #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.nh_0 > 0 and simulink_opts.inputs.uh_0 -%}  {#- uh_0 #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.nh_e > 0 and simulink_opts.inputs.lh_e -%}  {#- lh_e #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.nh_e > 0 and simulink_opts.inputs.uh_e -%}  {#- uh_e #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if dims.ny_0 > 0 and simulink_opts.inputs.cost_W_0 %}  {#- cost_W_0 #}
    {%- set n_inputs = n_inputs + 1 %}
  {%- endif -%}
  {%- if dims.ny > 0 and simulink_opts.inputs.cost_W %}  {#- cost_W #}
    {%- set n_inputs = n_inputs + 1 %}
  {%- endif -%}
  {%- if dims.ny_e > 0 and simulink_opts.inputs.cost_W_e %}  {#- cost_W_e #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}

  {%- if ns_total > 0 and simulink_opts.inputs.cost_zl %}  {#- cost_zl #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if ns_total > 0 and simulink_opts.inputs.cost_zu %}  {#- cost_zu #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if ns_total > 0 and simulink_opts.inputs.cost_Zl %}  {#- cost_Zl #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}
  {%- if ns_total > 0 and simulink_opts.inputs.cost_Zu %}  {#- cost_Zu #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}

  {%- if simulink_opts.inputs.reset_solver -%}  {#- reset_solver #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}

  {%- if simulink_opts.inputs.ignore_inits -%}  {#- ignore_inits #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}

  {%- if simulink_opts.inputs.x_init -%}  {#- x_init #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}

  {%- if simulink_opts.inputs.u_init -%}  {#- u_init #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}

  {%- if simulink_opts.inputs.pi_init -%}  {#- pi_init #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}

  {%- if simulink_opts.inputs.slacks_init -%}  {#- slacks_init #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}

  {%- if simulink_opts.inputs.rti_phase -%}  {#- rti_phase #}
    {%- set n_inputs = n_inputs + 1 -%}
  {%- endif -%}

{%- if simulink_opts.customizable_inputs %}
  {#- customizable inputs #}
  {%- for input_name, input_spec in simulink_opts.customizable_inputs -%}
    {%- if input_name is starting_with("sparse_parameter") -%}
      {%- set_global n_inputs = n_inputs + 1 -%}
    {%- else %}
      {{ throw(message = "only kind of supported customizable input are sparse_parameter, sparse_parameter_stagewise") }}
    {%- endif -%}
  {%- endfor -%}
{%- endif -%}

    // specify the number of input ports
    if ( !ssSetNumInputPorts(S, {{ n_inputs }}) )
        return;

    // specify the number of output ports
    {%- set_global n_outputs = 0 %}
    {%- for key, val in simulink_opts.outputs %}
      {%- if val == 1 %}
        {%- set_global n_outputs = n_outputs + val %}
      {%- elif val != 0 %}
        {{ throw(message = "simulink_opts.outputs must be 0 or 1, got val") }}
      {%- endif %}
    {%- endfor %}
    if ( !ssSetNumOutputPorts(S, {{ n_outputs }}) )
        return;

    // specify dimension information for the input ports
    {%- set i_input = -1 %}{# note here i_input is 0-based #}
  {%- if dims.nbx_0 > 0 and simulink_opts.inputs.lbx_0 -%}  {#- lbx_0 #}
    {%- set i_input = i_input + 1 %}
    // lbx_0
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nbx_0 }});
  {%- endif %}
  {%- if dims.nbx_0 > 0 and simulink_opts.inputs.ubx_0 -%}  {#- ubx_0 #}
    {%- set i_input = i_input + 1 %}
    // ubx_0
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nbx_0 }});
  {%- endif %}

  {%- if dims.np > 0 and simulink_opts.inputs.parameter_traj -%}  {#- parameter_traj #}
    {%- set i_input = i_input + 1 %}
    // parameters
    ssSetInputPortVectorDimension(S, {{ i_input }}, (N+1) * {{ dims.np }});
  {%- endif %}

  {%- if dims.ny > 0 and simulink_opts.inputs.y_ref_0 %}
    {%- set i_input = i_input + 1 %}
    // y_ref_0
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.ny_0 }});
  {%- endif %}

  {%- if dims.ny > 0 and solver_options.N_horizon > 1 and simulink_opts.inputs.y_ref %}
    {%- set i_input = i_input + 1 %}
    // y_ref
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ (solver_options.N_horizon-1) * dims.ny }});
  {%- endif %}

  {%- if dims.ny_e > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.y_ref_e %}
    {%- set i_input = i_input + 1 %}
    // y_ref_e
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.ny_e }});
  {%- endif %}

  {%- if dims.nbx > 0 and solver_options.N_horizon > 1 and simulink_opts.inputs.lbx -%}  {#- lbx #}
    {%- set i_input = i_input + 1 %}
    // lbx
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ (solver_options.N_horizon-1) * dims.nbx }});
  {%- endif %}
  {%- if dims.nbx > 0 and solver_options.N_horizon > 1 and simulink_opts.inputs.ubx -%}  {#- ubx #}
    {%- set i_input = i_input + 1 %}
    // ubx
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ (solver_options.N_horizon-1) * dims.nbx }});
  {%- endif %}

  {%- if dims.nbx_e > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.lbx_e -%}  {#- lbx_e #}
    {%- set i_input = i_input + 1 %}
    // lbx_e
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nbx_e }});
  {%- endif %}
  {%- if dims.nbx_e > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.ubx_e -%}  {#- ubx_e #}
    {%- set i_input = i_input + 1 %}
    // ubx_e
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nbx_e }});
  {%- endif %}

  {%- if dims.nbu > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.lbu -%}  {#- lbu #}
    {%- set i_input = i_input + 1 %}
    // lbu
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ solver_options.N_horizon*dims.nbu }});
  {%- endif -%}
  {%- if dims.nbu > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.ubu -%}  {#- ubu #}
    {%- set i_input = i_input + 1 %}
    // ubu
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ solver_options.N_horizon*dims.nbu }});
  {%- endif -%}


  {%- if dims.ng > 0 and simulink_opts.inputs.lg -%}  {#- lg #}
    {%- set i_input = i_input + 1 %}
    // lg
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ solver_options.N_horizon*dims.ng }});
  {%- endif -%}
  {%- if dims.ng > 0 and simulink_opts.inputs.ug -%}  {#- ug #}
    {%- set i_input = i_input + 1 %}
    // ug
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ solver_options.N_horizon*dims.ng }});
  {%- endif -%}

  {%- if dims.nh > 0 and simulink_opts.inputs.lh -%}  {#- lh #}
    {%- set i_input = i_input + 1 %}
    // lh
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ solver_options.N_horizon*dims.nh }});
  {%- endif -%}
  {%- if dims.nh > 0 and simulink_opts.inputs.uh -%}  {#- uh #}
    {%- set i_input = i_input + 1 %}
    // uh
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ solver_options.N_horizon*dims.nh }});
  {%- endif -%}

{%- if dims.nh_0 > 0 and simulink_opts.inputs.lh_0 -%}  {#- lh_0 #}
    {%- set i_input = i_input + 1 %}
    // lh_0
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nh_0 }});
  {%- endif -%}
  {%- if dims.nh_0 > 0 and simulink_opts.inputs.uh_0 -%}  {#- uh_0 #}
    {%- set i_input = i_input + 1 %}
    // uh_0
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nh_0 }});
  {%- endif -%}

  {%- if dims.nh_e > 0 and simulink_opts.inputs.lh_e -%}  {#- lh_e #}
    {%- set i_input = i_input + 1 %}
    // lh_e
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nh_e }});
  {%- endif -%}
  {%- if dims.nh_e > 0 and simulink_opts.inputs.uh_e -%}  {#- uh_e #}
    {%- set i_input = i_input + 1 %}
    // uh_e
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nh_e }});
  {%- endif -%}

  {%- if dims.ny_0 > 0 and simulink_opts.inputs.cost_W_0 %}  {#- cost_W_0 #}
    {%- set i_input = i_input + 1 %}
    // cost_W_0
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.ny_0 * dims.ny_0 }});
  {%- endif %}

  {%- if dims.ny > 0 and simulink_opts.inputs.cost_W %}  {#- cost_W #}
    {%- set i_input = i_input + 1 %}
    // cost_W
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.ny * dims.ny }});
  {%- endif %}

  {%- if dims.ny_e > 0 and simulink_opts.inputs.cost_W_e %}  {#- cost_W_e #}
    {%- set i_input = i_input + 1 %}
    // cost_W_e
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.ny_e * dims.ny_e }});
  {%- endif %}

  {%- if ns_total > 0 and simulink_opts.inputs.cost_zl %}  {#- cost_zl #}
    {%- set i_input = i_input + 1 %}
    // cost_zl
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ ns_total }});
  {%- endif %}
  {%- if ns_total > 0 and simulink_opts.inputs.cost_zu %}  {#- cost_zu #}
    {%- set i_input = i_input + 1 %}
    // cost_zu
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ ns_total }});
  {%- endif %}
  {%- if ns_total > 0 and simulink_opts.inputs.cost_Zl %}  {#- cost_Zl #}
    {%- set i_input = i_input + 1 %}
    // cost_Zl
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ ns_total }});
  {%- endif %}
  {%- if ns_total > 0 and simulink_opts.inputs.cost_Zu %}  {#- cost_Zu #}
    {%- set i_input = i_input + 1 %}
    // cost_Zu
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ ns_total }});
  {%- endif %}

  {%- if simulink_opts.inputs.reset_solver -%}  {#- reset_solver #}
    {%- set i_input = i_input + 1 %}
    // reset_solver
    ssSetInputPortVectorDimension(S, {{ i_input }}, 1);
  {%- endif -%}

  {%- if simulink_opts.inputs.ignore_inits -%}  {#- ignore_inits #}
    {%- set i_input = i_input + 1 %}
    // ignore_inits
    ssSetInputPortVectorDimension(S, {{ i_input }}, 1);
  {%- endif -%}

  {%- if simulink_opts.inputs.x_init -%}  {#- x_init #}
    {%- set i_input = i_input + 1 %}
    // x_init
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nx * (solver_options.N_horizon+1) }});
  {%- endif -%}

  {%- if simulink_opts.inputs.u_init -%}  {#- u_init #}
    {%- set i_input = i_input + 1 %}
    // u_init
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nu * (solver_options.N_horizon) }});
  {%- endif -%}

  {%- if simulink_opts.inputs.pi_init -%}  {#- pi_init #}
    {%- set i_input = i_input + 1 %}
    // pi_init
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ dims.nx * (solver_options.N_horizon) }});
  {%- endif -%}

  {%- if simulink_opts.inputs.slacks_init -%}  {#- slacks_init #}
    {%- set i_input = i_input + 1 %}
    // slacks_init
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ 2* ns_total }});
  {%- endif -%}


  {%- if simulink_opts.inputs.rti_phase -%}  {#- rti_phase #}
    {%- set i_input = i_input + 1 %}
    // rti_phase
    ssSetInputPortVectorDimension(S, {{ i_input }}, 1);
  {%- endif -%}


{%- if simulink_opts.customizable_inputs %}
  {#- customizable inputs #}
  {%- for input_name, input_spec in simulink_opts.customizable_inputs -%}
    {%- if input_name is starting_with("sparse_parameter") -%}
      {% set param_length = input_spec.parameter_indices | length %}
      {% set port_name = input_name | replace(from="sparse_parameter_", to="") %}
      {% set stage_idx_0 = input_spec.stage_idx_0 %}
      {% set stage_idx_e = input_spec.stage_idx_e %}
      {%- set_global i_input = i_input + 1 %}
    // {{ port_name }}
    ssSetInputPortVectorDimension(S, {{ i_input }}, {{ 1 + (stage_idx_e - stage_idx_0 + 1) * param_length }});
    {%- else %}
      {{ throw(message = "only kind of supported customizable input are sparse_parameter.") }}
    {%- endif -%}
  {%- endfor -%}
{%- endif -%}


    /* specify dimension information for the OUTPUT ports */
    {%- set i_output = -1 %}{# note here i_output is 0-based #}
  {%- if dims.nu > 0 and simulink_opts.outputs.u0 == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, {{ dims.nu }} );
  {%- endif %}

  {%- if simulink_opts.outputs.utraj == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, {{ dims.nu * solver_options.N_horizon }} );
  {%- endif %}

  {%- if simulink_opts.outputs.xtraj == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, {{ dims.nx * (solver_options.N_horizon+1) }} );
  {%- endif %}

  {%- if simulink_opts.outputs.ztraj == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, {{ dims.nz * solver_options.N_horizon }} );
  {%- endif %}

  {%- if simulink_opts.outputs.pi_all == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, {{ dims.nx * solver_options.N_horizon }} );
  {%- endif %}

  {%- if simulink_opts.outputs.slack_values == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, 2*{{ ns_total }} );
  {%- endif %}

  {%- if simulink_opts.outputs.solver_status == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, 1 );
  {%- endif %}

  {%- if simulink_opts.outputs.cost_value == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, 1 );
  {%- endif %}

  {%- if simulink_opts.outputs.KKT_residual == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, 1 );
  {%- endif %}

  {%- if simulink_opts.outputs.KKT_residuals == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, 4 );
  {%- endif %}

  {%- if solver_options.N_horizon > 0 and simulink_opts.outputs.x1 == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, {{ dims.nx }} ); // state at shooting node 1
  {%- endif %}

  {%- if simulink_opts.outputs.CPU_time == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, 1);
  {%- endif %}

  {%- if simulink_opts.outputs.CPU_time_sim == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, 1);
  {%- endif %}

  {%- if simulink_opts.outputs.CPU_time_qp == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, 1);
  {%- endif %}

  {%- if simulink_opts.outputs.CPU_time_lin == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, 1);
  {%- endif %}

  {%- if simulink_opts.outputs.sqp_iter == 1 %}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, 1 );
  {%- endif %}
  {%- if simulink_opts.outputs.parameter_traj -%}  {#- parameter_traj #}
    {%- set i_output = i_output + 1 %}
    ssSetOutputPortVectorDimension(S, {{ i_output }}, {{ dims.np * (solver_options.N_horizon + 1) }});
  {%- endif -%}

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
    {{ name }}_solver_capsule *capsule = {{ name }}_acados_create_capsule();
    {{ name }}_acados_create(capsule);

    ssSetUserData(S, (void*)capsule);
}


static void mdlOutputs(SimStruct *S, int_T tid)
{
    {{ name }}_solver_capsule *capsule = ssGetUserData(S);
    ocp_nlp_config *nlp_config = {{ name }}_acados_get_nlp_config(capsule);
    ocp_nlp_dims *nlp_dims = {{ name }}_acados_get_nlp_dims(capsule);
    ocp_nlp_in *nlp_in = {{ name }}_acados_get_nlp_in(capsule);
    ocp_nlp_out *nlp_out = {{ name }}_acados_get_nlp_out(capsule);

    InputRealPtrsType in_sign;

    int N = {{ model.name | upper }}_N;

    {%- set buffer_sizes = [dims.nx, dims.nu, dims.nbx_0, dims.np, dims.nbx, dims.nbx_e, dims.nbu, dims.ng, dims.nh, dims.nh_0, dims.ng_e, dims.nh_e, dims.ns, dims.ns_0, dims.ns_e] -%}

  {%- if dims.ny_0 > 0 and simulink_opts.inputs.y_ref_0 %}  {# y_ref_0 #}
    {%- set buffer_sizes = buffer_sizes | concat(with=(dims.ny_0)) %}
  {%- endif %}
  {%- if dims.ny > 0 and solver_options.N_horizon > 1 and simulink_opts.inputs.y_ref %}  {# y_ref #}
    {%- set buffer_sizes = buffer_sizes | concat(with=(dims.ny)) %}
  {%- endif %}
  {%- if dims.ny_e > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.y_ref_e %}  {# y_ref_e #}
    {%- set buffer_sizes = buffer_sizes | concat(with=(dims.ny_e)) %}
  {%- endif %}

  {%- if dims.ny_0 > 0 and simulink_opts.inputs.cost_W_0 %}  {#- cost_W_0 #}
    {%- set buffer_sizes = buffer_sizes | concat(with=(dims.ny_0 * dims.ny_0)) %}
  {%- endif %}
  {%- if dims.ny > 0 and simulink_opts.inputs.cost_W %}  {#- cost_W #}
    {%- set buffer_sizes = buffer_sizes | concat(with=(dims.ny * dims.ny)) %}
  {%- endif %}
  {%- if dims.ny_e > 0 and simulink_opts.inputs.cost_W_e %}  {#- cost_W_e #}
    {%- set buffer_sizes = buffer_sizes | concat(with=(dims.ny_e * dims.ny_e)) %}
  {%- endif %}

    // local buffer
    {%- set buffer_size = buffer_sizes | sort | last %}
    double buffer[{{ buffer_size }}];
    double tmp_double;
    int tmp_offset, tmp_int;

    /* go through inputs */
    {%- set i_input = -1 %}
  {%- if dims.nbx_0 > 0 and simulink_opts.inputs.lbx_0 -%}  {#- lbx_0 #}
    // lbx_0
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int i = 0; i < {{ dims.nbx_0 }}; i++)
        buffer[i] = (double)(*in_sign[i]);

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", buffer);
  {%- endif %}

  {%- if dims.nbx_0 > 0 and simulink_opts.inputs.ubx_0 -%}  {#- ubx_0 #}
    // ubx_0
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int i = 0; i < {{ dims.nbx_0 }}; i++)
        buffer[i] = (double)(*in_sign[i]);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", buffer);
  {%- endif %}

  {%- if dims.np > 0 and simulink_opts.inputs.parameter_traj -%}  {#- parameter_traj #}
    // parameter_traj
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    // update value of parameters
    tmp_offset = 0;
    for (int stage = 0; stage <= N; stage++)
    {
        tmp_int = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, stage, "p");
        for (int jj = 0; jj < tmp_int; jj++)
        {
            buffer[jj] = (double)(*in_sign[tmp_offset+jj]);
        }
        {{ name }}_acados_update_params(capsule, stage, buffer, tmp_int);
        tmp_offset += tmp_int;
    }
  {%- endif %}

  {% if dims.ny_0 > 0 and simulink_opts.inputs.y_ref_0 %}
    // y_ref_0
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});

    for (int i = 0; i < {{ dims.ny_0 }}; i++)
        buffer[i] = (double)(*in_sign[i]);

    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, 0, "yref", (void *) buffer);
  {%- endif %}

  {% if dims.ny > 0 and solver_options.N_horizon > 1 and simulink_opts.inputs.y_ref %}
    // y_ref - for stages 1 to N-1
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});

    for (int stage = 1; stage < N; stage++)
    {
        for (int jj = 0; jj < {{ dims.ny }}; jj++)
            buffer[jj] = (double)(*in_sign[(stage-1)*{{ dims.ny }}+jj]);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, stage, "yref", (void *) buffer);
    }
  {%- endif %}

  {% if dims.ny_e > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.y_ref_e %}
    // y_ref_e
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});

    for (int i = 0; i < {{ dims.ny_e }}; i++)
        buffer[i] = (double)(*in_sign[i]);

    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "yref", (void *) buffer);
  {%- endif %}

  {%- if dims.nbx > 0 and solver_options.N_horizon > 1 and simulink_opts.inputs.lbx -%}  {#- lbx #}
    // lbx
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int stage = 1; stage < N; stage++)
    {
        for (int jj = 0; jj < {{ dims.nbx }}; jj++)
            buffer[jj] = (double)(*in_sign[(stage-1)*{{ dims.nbx }}+jj]);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, stage, "lbx", (void *) buffer);
    }
  {%- endif %}
  {%- if dims.nbx > 0 and solver_options.N_horizon > 1 and simulink_opts.inputs.ubx -%}  {#- ubx #}
    // ubx
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int stage = 1; stage < N; stage++)
    {
        for (int jj = 0; jj < {{ dims.nbx }}; jj++)
            buffer[jj] = (double)(*in_sign[(stage-1)*{{ dims.nbx }}+jj]);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, stage, "ubx", (void *) buffer);
    }
  {%- endif %}


  {%- if dims.nbx_e > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.lbx_e -%}  {#- lbx_e #}
    // lbx_e
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});

    for (int i = 0; i < {{ dims.nbx_e }}; i++)
        buffer[i] = (double)(*in_sign[i]);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lbx", buffer);
  {%- endif %}
  {%- if dims.nbx_e > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.ubx_e -%}  {#- ubx_e #}
    // ubx_e
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});

    for (int i = 0; i < {{ dims.nbx_e }}; i++)
        buffer[i] = (double)(*in_sign[i]);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "ubx", buffer);
  {%- endif %}


  {%- if dims.nbu > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.lbu -%}  {#- lbu #}
    // lbu
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int stage = 0; stage < N; stage++)
    {
        for (int jj = 0; jj < {{ dims.nbu }}; jj++)
            buffer[jj] = (double)(*in_sign[stage*{{ dims.nbu }}+jj]);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, stage, "lbu", (void *) buffer);
    }
  {%- endif -%}
  {%- if dims.nbu > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.ubu -%}  {#- ubu #}
    // ubu
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int stage = 0; stage < N; stage++)
    {
        for (int jj = 0; jj < {{ dims.nbu }}; jj++)
            buffer[jj] = (double)(*in_sign[stage*{{ dims.nbu }}+jj]);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, stage, "ubu", (void *) buffer);
    }
  {%- endif -%}

  {%- if dims.ng > 0 and simulink_opts.inputs.lg -%}  {#- lg #}
    // lg
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});

    for (int stage = 0; stage < N; stage++)
    {
        for (int jj = 0; jj < {{ dims.ng }}; jj++)
            buffer[jj] = (double)(*in_sign[stage*{{ dims.ng }}+jj]);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, stage, "lg", (void *) buffer);
    }
  {%- endif -%}
  {%- if dims.ng > 0 and simulink_opts.inputs.ug -%}  {#- ug #}
    // ug
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});

    for (int stage = 0; stage < N; stage++)
    {
        for (int jj = 0; jj < {{ dims.ng }}; jj++)
            buffer[jj] = (double)(*in_sign[stage*{{ dims.ng }}+jj]);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, stage, "ug", (void *) buffer);
    }
  {%- endif -%}

  {%- if dims.nh > 0 and simulink_opts.inputs.lh -%}  {#- lh #}
    // lh
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});

    for (int stage = 0; stage < N; stage++)
    {
        for (int jj = 0; jj < {{ dims.nh }}; jj++)
            buffer[jj] = (double)(*in_sign[stage*{{ dims.nh }}+jj]);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, stage, "lh", (void *) buffer);
    }
  {%- endif -%}
  {%- if dims.nh > 0 and simulink_opts.inputs.uh -%}  {#- uh #}
    // uh
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});

    for (int stage = 0; stage < N; stage++)
    {
        for (int jj = 0; jj < {{ dims.nh }}; jj++)
            buffer[jj] = (double)(*in_sign[stage*{{ dims.nh }}+jj]);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, stage, "uh", (void *) buffer);
    }
  {%- endif -%}

{%- if dims.nh_0 > 0 and simulink_opts.inputs.lh_0 -%}  {#- lh_0 #}
    // lh_0
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int i = 0; i < {{ dims.nh_0 }}; i++)
        buffer[i] = (double)(*in_sign[i]);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lh", buffer);
  {%- endif -%}
  {%- if dims.nh_0 > 0 and simulink_opts.inputs.uh_0 -%}  {#- uh_0 #}
    // uh_0
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int i = 0; i < {{ dims.nh_0 }}; i++)
        buffer[i] = (double)(*in_sign[i]);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "uh", buffer);
  {%- endif -%}

  {%- if dims.nh_e > 0 and simulink_opts.inputs.lh_e -%}  {#- lh_e #}
    // lh_e
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int i = 0; i < {{ dims.nh_e }}; i++)
        buffer[i] = (double)(*in_sign[i]);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lh", buffer);
  {%- endif -%}
  {%- if dims.nh_e > 0 and simulink_opts.inputs.uh_e -%}  {#- uh_e #}
    // uh_e
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int i = 0; i < {{ dims.nh_e }}; i++)
        buffer[i] = (double)(*in_sign[i]);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "uh", buffer);
  {%- endif -%}

  {%- if dims.ny_0 > 0 and simulink_opts.inputs.cost_W_0 %}  {#- cost_W_0 #}
    // cost_W_0
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int i = 0; i < {{ dims.ny_0 * dims.ny_0 }}; i++)
        buffer[i] = (double)(*in_sign[i]);

    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, 0, "W", buffer);
  {%- endif %}

  {%- if dims.ny > 0 and simulink_opts.inputs.cost_W %}  {#- cost_W #}
    // cost_W
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int i = 0; i < {{ dims.ny * dims.ny }}; i++)
        buffer[i] = (double)(*in_sign[i]);

    for (int stage = 1; stage < N; stage++)
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, stage, "W", buffer);
  {%- endif %}

  {%- if dims.ny_e > 0 and simulink_opts.inputs.cost_W_e %}  {#- cost_W_e #}
    // cost_W_e
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    for (int i = 0; i < {{ dims.ny_e * dims.ny_e }}; i++)
        buffer[i] = (double)(*in_sign[i]);

    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "W", buffer);
  {%- endif %}

  {%- if ns_total > 0 and simulink_opts.inputs.cost_zl %}  {#- cost_zl #}
    // cost_zl
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    tmp_offset = 0;
    for (int stage = 0; stage <= N; stage++)
    {
        tmp_int = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, stage, "zl");
        for (int i = 0; i < tmp_int; i++)
            buffer[i] = (double)(*in_sign[tmp_offset+i]);
        tmp_offset += tmp_int;
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, stage, "zl", buffer);
    }
  {%- endif %}
  {%- if ns_total > 0 and simulink_opts.inputs.cost_zl %}  {#- cost_zu #}
    // cost_zu
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    tmp_offset = 0;
    for (int stage = 0; stage <= N; stage++)
    {
        tmp_int = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, stage, "zu");
        for (int i = 0; i < tmp_int; i++)
            buffer[i] = (double)(*in_sign[tmp_offset+i]);
        tmp_offset += tmp_int;
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, stage, "zu", buffer);
    }
  {%- endif %}

  {%- if ns_total > 0 and simulink_opts.inputs.cost_Zl %}  {#- cost_Zl #}
    // cost_Zl
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    tmp_offset = 0;
    for (int stage = 0; stage <= N; stage++)
    {
        tmp_int = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, stage, "Zl");
        for (int i = 0; i < tmp_int; i++)
            buffer[i] = (double)(*in_sign[tmp_offset+i]);
        tmp_offset += tmp_int;
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, stage, "Zl", buffer);
    }
  {%- endif %}
  {%- if ns_total > 0 and simulink_opts.inputs.cost_Zl %}  {#- cost_Zu #}
    // cost_Zu
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    tmp_offset = 0;
    for (int stage = 0; stage <= N; stage++)
    {
        tmp_int = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, stage, "Zu");
        for (int i = 0; i < tmp_int; i++)
            buffer[i] = (double)(*in_sign[tmp_offset+i]);
        tmp_offset += tmp_int;
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, stage, "Zu", buffer);
    }
  {%- endif %}

  {%- if simulink_opts.inputs.reset_solver %}  {#- reset_solver #}
    // reset_solver
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    double reset = (double)(*in_sign[0]);
    if (reset)
    {
        {{ name }}_acados_reset(capsule, 1);
    }
  {%- endif %}

    int ignore_inits = 0;
  {%- if simulink_opts.inputs.ignore_inits %}  {#- ignore_inits #}
    // ignore_inits
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    ignore_inits = (int)(*in_sign[0]);
  {%- endif %}
    // ssPrintf("ignore_inits = %d\n", ignore_inits);

    if (ignore_inits == 0)
    {
      {%- if simulink_opts.inputs.x_init %}  {#- x_init #}
        // x_init
        {%- set i_input = i_input + 1 %}
        in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
        tmp_offset = 0;
        for (int stage = 0; stage < {{ solver_options.N_horizon + 1 }}; stage++)
        {
            tmp_int = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, stage, "x");
            for (int jj = 0; jj < tmp_int; jj++)
                buffer[jj] = (double)(*in_sign[tmp_offset+jj]);
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, stage, "x", (void *) buffer);
            tmp_offset += tmp_int;
        }
      {%- endif %}

      {%- if simulink_opts.inputs.u_init %}  {#- u_init #}
        // u_init
        {%- set i_input = i_input + 1 %}
        in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
        tmp_offset = 0;
        for (int stage = 0; stage < N; stage++)
        {
            tmp_int = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, stage, "u");
            for (int jj = 0; jj < {{ dims.nu }}; jj++)
                buffer[jj] = (double)(*in_sign[tmp_offset+jj]);
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, stage, "u", (void *) buffer);
            tmp_offset += tmp_int;
        }
      {%- endif %}

      {%- if simulink_opts.inputs.pi_init %}  {#- pi_init #}
        // pi_init
        {%- set i_input = i_input + 1 %}
        in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
        tmp_offset = 0;
        for (int stage = 0; stage < N; stage++)
        {
            tmp_int = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, stage, "pi");
            for (int jj = 0; jj < tmp_int; jj++)
                buffer[jj] = (double)(*in_sign[tmp_offset+jj]);
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, stage, "pi", (void *) buffer);
            tmp_offset += tmp_int;
        }
      {%- endif %}

      {%- if simulink_opts.inputs.slacks_init %}  {#- slacks_init #}
        // slacks_init
        {%- set i_input = i_input + 1 %}
        in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
        tmp_offset = 0;
        for (int stage = 0; stage <= N; stage++)
        {
            tmp_int = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, stage, "sl");
            // set sl
            for (int jj = 0; jj < tmp_int; jj++)
                buffer[jj] = (double)(*in_sign[tmp_offset+jj]);
            tmp_offset += tmp_int;
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, stage, "sl", (void *) buffer);
            // set su
            for (int jj = 0; jj < tmp_int; jj++)
                buffer[jj] = (double)(*in_sign[tmp_offset+jj]);
            tmp_offset += tmp_int;
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, stage, "su", (void *) buffer);
        }
      {%- endif %}
    }

  {%- if simulink_opts.inputs.rti_phase %}  {#- rti_phase #}
    {%- set i_input = i_input + 1 %}
    in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});
    double rti_phase_double = (double)(*in_sign[0]);
    int rti_phase = (int) rti_phase_double;

    ocp_nlp_solver_opts_set(nlp_config, capsule->nlp_opts, "rti_phase", &rti_phase);
  {%- endif %}


{%- if simulink_opts.customizable_inputs %}
  {#- customizable inputs #}
  {%- for input_name, input_spec in simulink_opts.customizable_inputs -%}
    {%- if input_name is starting_with("sparse_parameter") %}
    // length of parameter_indices {{ input_spec.parameter_indices | length }}
    {%- set_global i_input = i_input + 1 %}
    {% set param_length = input_spec.parameter_indices | length %}
    {% set port_name = input_name | replace(from="sparse_parameter_", to="") %}
    {% set stage_idx_0 = input_spec.stage_idx_0 %}
    {% set stage_idx_e = input_spec.stage_idx_e %}
    // {{ port_name }}

        in_sign = ssGetInputPortRealSignalPtrs(S, {{ i_input }});

        tmp_double = (double)(*in_sign[0]); // decides if update is done.
        if (tmp_double)
        {
            int idx[{{ param_length }}];
            {% for item in input_spec.parameter_indices %}
            idx[{{ loop.index0 }}] = {{ item }};
            {%- endfor %}

            // update for stages
            for (int stage = {{ stage_idx_0 }}; stage < {{ stage_idx_e }}+1; stage++)
            {
                tmp_offset = 1 + (stage - {{ stage_idx_0 }}) * {{ param_length }};
                // copy new parameter values to buffer
                for (int jj = 0; jj < {{ param_length }}; jj++)
                {
                    buffer[jj] = (double)(*in_sign[jj + tmp_offset]);
                }
                {{ name }}_acados_update_params_sparse(capsule, stage, idx, buffer, {{ param_length }});
            }
        }
    {%- endif -%}
  {%- endfor -%}
{%- endif -%}

    /* call solver */
  {%- if custom_update_filename == "" and not simulink_opts.inputs.rti_phase %}
    int acados_status = {{ name }}_acados_solve(capsule);
    // get time
    ocp_nlp_get(nlp_config, capsule->nlp_solver, "time_tot", (void *) buffer);
    tmp_double = buffer[0];

  {%- elif simulink_opts.inputs.rti_phase %}{# SPLIT RTI PHASE#}
    {% if solver_options.nlp_solver_type != "SQP_RTI" %}
    rti_phase input only supported for nlp_solver_type == "SQP_RTI"!
    {% elif custom_update_filename != "" %}
    rti_phase input only supported for custom_update_filename == ""!
    {% else %}
    ocp_nlp_solver_opts_set(nlp_config, capsule->nlp_opts, "rti_phase", &rti_phase);
    int acados_status = {{ name }}_acados_solve(capsule);
    // get time
    ocp_nlp_get(nlp_config, capsule->nlp_solver, "time_tot", (void *) buffer);
    tmp_double = buffer[0];
    {%- endif %}
  {%- elif solver_options.nlp_solver_type == "SQP_RTI" %}{# if custom_update_filename != "" #}
    // preparation
    int rti_phase = 1;
    ocp_nlp_solver_opts_set(nlp_config, capsule->nlp_opts, "rti_phase", &rti_phase);
    int acados_status = {{ name }}_acados_solve(capsule);

    // preparation time
    ocp_nlp_get(nlp_config, capsule->nlp_solver, "time_tot", (void *) buffer);
    tmp_double = buffer[0];

    // call custom update function
    int data_len = 0;
    double* c_data; // TODO: only works with empty..
    acados_status = {{ name }}_acados_custom_update(capsule, c_data, data_len);

    // feedback
    rti_phase = 2;
    ocp_nlp_solver_opts_set(nlp_config, capsule->nlp_opts, "rti_phase", &rti_phase);
    acados_status = {{ name }}_acados_solve(capsule);
    // feedback time
    ocp_nlp_get(nlp_config, capsule->nlp_solver, "time_tot", (void *) buffer);
    tmp_double += buffer[0];
  {%- else -%}
    Simulink block with custom solver template only works with SQP_RTI!
  {%- endif %}

    /* set outputs */
    double *out_ptr;
    {%- set i_output = -1 -%}{# note here i_output is 0-based #}
  {%- if dims.nu > 0 and simulink_opts.outputs.u0 == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 0, "u", (void *) out_ptr);
  {%- endif %}

  {%- if simulink_opts.outputs.utraj == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    tmp_offset = 0;
    for (int stage = 0; stage < N; stage++)
    {
        tmp_int = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, stage, "u");
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, stage, "u", (void *) (out_ptr + tmp_offset));
        tmp_offset += tmp_int;
    }
  {%- endif %}

  {% if simulink_opts.outputs.xtraj == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    tmp_offset = 0;
    for (int stage = 0; stage < {{ solver_options.N_horizon + 1 }}; stage++)
    {
        tmp_int = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, stage, "x");
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, stage, "x", (void *) (out_ptr + tmp_offset));
        tmp_offset += tmp_int;
    }
  {%- endif %}

  {% if simulink_opts.outputs.ztraj == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    tmp_offset = 0;
    for (int stage = 0; stage < N; stage++)
    {
        tmp_int = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, stage, "z");
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, stage, "z", (void *) (out_ptr + tmp_offset));
        tmp_offset += tmp_int;
    }
  {%- endif %}

  {% if simulink_opts.outputs.pi_all == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    tmp_offset = 0;
    for (int stage = 0; stage < N; stage++)
    {
        tmp_int = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, stage, "pi");
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, stage, "pi", (void *) (out_ptr + stage * {{ dims.nx }}));
        tmp_offset += tmp_int;
    }
  {%- endif %}

  {% if simulink_opts.outputs.slack_values == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    tmp_offset = 0;
    for (int stage = 0; stage <= N; stage++)
    {
        tmp_int = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, stage, "sl");
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, stage, "sl", (void *) (out_ptr + tmp_offset));
        tmp_offset += tmp_int;
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, stage, "su", (void *) (out_ptr + tmp_offset));
        tmp_offset += tmp_int;
    }
  {%- endif %}

  {%- if simulink_opts.outputs.solver_status == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    *out_ptr = (double) acados_status;
  {%- endif %}

  {%- if simulink_opts.outputs.cost_value == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    ocp_nlp_eval_cost(capsule->nlp_solver, nlp_in, nlp_out);
    ocp_nlp_get(nlp_config, capsule->nlp_solver, "cost_value", (void *) out_ptr);
  {%- endif %}

  {%- if simulink_opts.outputs.KKT_residual == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    *out_ptr = (double) nlp_out->inf_norm_res;
  {%- endif %}

  {%- if simulink_opts.outputs.KKT_residuals == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});

    {%- if solver_options.nlp_solver_type == "SQP_RTI" %}
    ocp_nlp_eval_residuals(capsule->nlp_solver, nlp_in, nlp_out);
    {%- endif %}
    ocp_nlp_get(nlp_config, capsule->nlp_solver, "res_stat", (void *) &out_ptr[0]);
    ocp_nlp_get(nlp_config, capsule->nlp_solver, "res_eq", (void *) &out_ptr[1]);
    ocp_nlp_get(nlp_config, capsule->nlp_solver, "res_ineq", (void *) &out_ptr[2]);
    ocp_nlp_get(nlp_config, capsule->nlp_solver, "res_comp", (void *) &out_ptr[3]);
  {%- endif %}

  {%- if solver_options.N_horizon > 0 and simulink_opts.outputs.x1 == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 1, "x", (void *) out_ptr);
  {%- endif %}

  {%- if simulink_opts.outputs.CPU_time == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    out_ptr[0] = tmp_double;
  {%- endif -%}

  {%- if simulink_opts.outputs.CPU_time_sim == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    ocp_nlp_get(nlp_config, capsule->nlp_solver, "time_sim", (void *) out_ptr);
  {%- endif -%}

  {%- if simulink_opts.outputs.CPU_time_qp == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    ocp_nlp_get(nlp_config, capsule->nlp_solver, "time_qp", (void *) out_ptr);
  {%- endif -%}

  {%- if simulink_opts.outputs.CPU_time_lin == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    ocp_nlp_get(nlp_config, capsule->nlp_solver, "time_lin", (void *) out_ptr);
  {%- endif -%}

  {%- if simulink_opts.outputs.sqp_iter == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    // get sqp iter
    ocp_nlp_get(nlp_config, capsule->nlp_solver, "sqp_iter", (void *) &tmp_int);
    *out_ptr = (double) tmp_int;
  {%- endif %}

  {% if simulink_opts.outputs.parameter_traj == 1 %}
    {%- set i_output = i_output + 1 %}
    out_ptr = ssGetOutputPortRealSignal(S, {{ i_output }});
    tmp_offset = 0;
    for (int stage = 0; stage < N+1; stage++)
    {
        tmp_int = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, stage, "p");
        ocp_nlp_in_get(nlp_config, nlp_dims, nlp_in, stage, "p", (void *) (out_ptr + tmp_offset));
        tmp_offset += tmp_int;
    }
  {%- endif %}

}

static void mdlTerminate(SimStruct *S)
{
    {{ name }}_solver_capsule *capsule = ssGetUserData(S);

    {{ name }}_acados_free(capsule);
    {{ name }}_acados_free_capsule(capsule);
}


#ifdef  MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
