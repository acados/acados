%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;

%
{%- if not solver_options.custom_update_filename %}
    {%- set custom_update_filename = "" %}
{% else %}
    {%- set custom_update_filename = solver_options.custom_update_filename %}
{%- endif %}


{% if problem_class == "OCP" %}
  {% set dims_e = dims %}
  {% set dims_0 = dims %}

  {%- set ns_total = dims.ns_0 + dims.ns_e + (solver_options.N_horizon - 1) * dims.ns %}
  {% set_global nx_total = dims.nx * (solver_options.N_horizon+1) %}
  {% set nbx_total = dims.nbx * (solver_options.N_horizon-1) %}{# Note: without initial and terminal node #}
  {% set nh_total = dims.nh * (solver_options.N_horizon-1) %}{# Note: without initial and terminal node #}
  {% set_global nu_total = dims.nu * (solver_options.N_horizon) %}
  {% set_global nbu_total = dims.nbu * (solver_options.N_horizon) %}
  {% set np_total = dims.np * (solver_options.N_horizon+1) %}
  {% set npi_total = dims.nx * (solver_options.N_horizon) %}
  {% set np_max = dims.np %}
  {% set nx_max = dims.nx %}
  {% set nu_max = dims.nu %}
  {% set ns_values = [dims.ns_0, dims.ns, dims.ns_e] %}
  {%- set ns_max = ns_values | sort | last %}
{% else %}
  {% set dims_0 = phases_dims | first %}
  {% set cost_0 = cost | first %}
  {% set constraints_0 = constraints | first %}
  {% set model_0 = model | first %}

  {% set cost_e = cost | last %}
  {% set constraints_e = constraints | last %}
  {% set dims_e = phases_dims | last %}
  {% set model_e = model | last %}

  {% set ns_total = dims_0.ns_0 %}
  {% set nx_total = 0 %}
  {% set nbx_total = 0 %}{# Note: without initial and terminal node #}
  {% set nh_total = 0 %}{# Note: without initial and terminal node #}
  {% set nu_total = 0 %}
  {% set nbu_total = 0 %}
  {% set nz_total = 0 %}
  {% set np_total = 0 %}
  {% set npi_total = 0 %}
  {% for jj in range(end=n_phases) %}{# phases loop !#}
    {% set_global ns_total = ns_total + (end_idx[jj] - cost_start_idx[jj]) * phases_dims[jj].ns %}
    {% set_global nx_total = nx_total + (end_idx[jj] - start_idx[jj]) * phases_dims[jj].nx %}
    {% set_global nbx_total = nbx_total + (end_idx[jj] - cost_start_idx[jj]) * phases_dims[jj].nbx %}
    {% set_global nh_total = nh_total + (end_idx[jj] - cost_start_idx[jj]) * phases_dims[jj].nh %}
    {% set_global nu_total = nu_total + (end_idx[jj] - start_idx[jj]) * phases_dims[jj].nu %}
    {% set_global nbu_total = nbu_total + (end_idx[jj] - start_idx[jj]) * phases_dims[jj].nbu %}
    {% set_global nz_total = nz_total + (end_idx[jj] - start_idx[jj]) * phases_dims[jj].nz %}
    {% set_global np_total = np_total + (end_idx[jj] - start_idx[jj]) * phases_dims[jj].np %}
    {% set_global npi_total = npi_total + (end_idx[jj] - start_idx[jj]) * phases_dims[jj].nx_next %}
  {% endfor %}{# phases loop !#}

  {% set_global nx_total = nx_total + dims_e.nx %}
  {% set_global np_total = np_total + dims_e.np %}
  {% set_global ns_total = ns_total + dims_e.ns_e %}

  {%- set nx_values = [] -%}
  {%- for jj in range(end=n_phases) %}
      {%- set_global nx_values = nx_values | concat(with=(phases_dims[jj].nx)) %}
  {%- endfor %}
  {%- set nx_max = nx_values | sort | last %}

  {%- set nu_values = [] -%}
  {%- for jj in range(end=n_phases) %}
      {%- set_global nu_values = nu_values | concat(with=(phases_dims[jj].nu)) %}
  {%- endfor %}
  {%- set nu_max = nu_values | sort | last %}

  {%- set np_values = [] -%}
  {%- for jj in range(end=n_phases) %}
      {%- set_global np_values = np_values | concat(with=(phases_dims[jj].np)) %}
  {%- endfor %}
  {%- set np_max = np_values | sort | last %}

  {%- set ns_values = [] -%}
  {%- for jj in range(end=n_phases) %}
      {%- set_global ns_values = ns_values | concat(with=(phases_dims[jj].ns)) %}
      {%- set_global ns_values = ns_values | concat(with=(phases_dims[jj].ns_0)) %}
      {%- set_global ns_values = ns_values | concat(with=(phases_dims[jj].ns_e)) %}
  {%- endfor %}
  {%- set ns_max = ns_values | sort | last %}
{%- endif %}

{# two brackets in math expression are not allowed by currently used tera #}
{%- set two_ns_total = 2 * ns_total %}

SOURCES = { ...
{%- for filename in external_function_files_model %}
            '{{ filename }}', ...
{%- endfor %}
{%- for filename in external_function_files_ocp %}
            '{{ filename }}', ...
{%- endfor %}
        {%- if custom_update_filename != "" %}
            '{{ custom_update_filename }}', ...
        {%- endif %}
            'acados_solver_sfunction_{{ name }}.c', ...
            'acados_solver_{{ name }}.c'
          };

INC_PATH = '{{ acados_include_path }}';

INCS = {['-I', fullfile(INC_PATH, 'blasfeo', 'include')], ...
        ['-I', fullfile(INC_PATH, 'hpipm', 'include')], ...
        ['-I', fullfile(INC_PATH, 'acados')], ...
        ['-I', fullfile(INC_PATH)]};

{% if solver_options.qp_solver is containing("QPOASES") %}
INCS{end+1} = ['-I', fullfile(INC_PATH, 'qpOASES_e')];
{% endif %}

CFLAGS = 'CFLAGS=$CFLAGS';
LDFLAGS = 'LDFLAGS=$LDFLAGS';
COMPFLAGS = 'COMPFLAGS=$COMPFLAGS';
COMPDEFINES = 'COMPDEFINES=$COMPDEFINES';

{% if solver_options.qp_solver is containing("QPOASES") %}
CFLAGS = [ CFLAGS, ' -DACADOS_WITH_QPOASES ' ];
COMPDEFINES = [ COMPDEFINES, ' -DACADOS_WITH_QPOASES ' ];
{%- elif solver_options.qp_solver is containing("OSQP") %}
CFLAGS = [ CFLAGS, ' -DACADOS_WITH_OSQP ' ];
COMPDEFINES = [ COMPDEFINES, ' -DACADOS_WITH_OSQP ' ];
{%- elif solver_options.qp_solver is containing("QPDUNES") %}
CFLAGS = [ CFLAGS, ' -DACADOS_WITH_QPDUNES ' ];
COMPDEFINES = [ COMPDEFINES, ' -DACADOS_WITH_QPDUNES ' ];
{%- elif solver_options.qp_solver is containing("DAQP") %}
CFLAGS = [ CFLAGS, ' -DACADOS_WITH_DAQP' ];
COMPDEFINES = [ COMPDEFINES, ' -DACADOS_WITH_DAQP' ];
{%- elif solver_options.qp_solver is containing("HPMPC") %}
CFLAGS = [ CFLAGS, ' -DACADOS_WITH_HPMPC ' ];
COMPDEFINES = [ COMPDEFINES, ' -DACADOS_WITH_HPMPC ' ];
{% endif %}

LIB_PATH = ['-L', fullfile('{{ acados_lib_path }}')];

LIBS = {'-lacados', '-lhpipm', '-lblasfeo'};

% acados linking libraries and flags
{%- if acados_link_libs and os and os == "pc" %}
LDFLAGS = [LDFLAGS ' {{ acados_link_libs.openmp }}'];
COMPFLAGS = [COMPFLAGS ' {{ acados_link_libs.openmp }}'];
LIBS{end+1} = '{{ acados_link_libs.qpoases }}';
LIBS{end+1} = '{{ acados_link_libs.hpmpc }}';
LIBS{end+1} = '{{ acados_link_libs.osqp }}';
{%- else %}
    {% if solver_options.qp_solver is containing("QPOASES") %}
LIBS{end+1} = '-lqpOASES_e';
    {% endif %}
    {% if solver_options.qp_solver is containing("DAQP") %}
LIBS{end+1} = '-ldaqp';
    {% endif %}
{%- endif %}

COMPFLAGS = [COMPFLAGS ' {{ solver_options.ext_fun_compile_flags }}'];
CFLAGS = [CFLAGS ' {{ solver_options.ext_fun_compile_flags }}'];

try
    %     mex('-v', '-O', CFLAGS, LDFLAGS, COMPFLAGS, COMPDEFINES, INCS{:}, ...
    mex('-O', CFLAGS, LDFLAGS, COMPFLAGS, COMPDEFINES, INCS{:}, ...
            LIB_PATH, LIBS{:}, SOURCES{:}, ...
            '-output', 'acados_solver_sfunction_{{ name }}' );
catch exception
    disp('make_sfun failed with the following exception:')
    disp(exception);
    disp(exception.message );
    disp('Try adding -v to the mex command above to get more information.')
    keyboard
end

fprintf( [ '\n\nSuccessfully created sfunction:\nacados_solver_sfunction_{{ name }}', '.', ...
    eval('mexext')] );


%% print note on usage of s-function, and create I/O port names vectors
fprintf('\n\nNote: Usage of Sfunction is as follows:\n')
input_note = 'Inputs are:\n';
i_in = 1;

global sfun_input_names
sfun_input_names = {};

{%- if dims_0.nbx_0 > 0 and simulink_opts.inputs.lbx_0 -%}  {#- lbx_0 #}
input_note = strcat(input_note, num2str(i_in), ') lbx_0 - lower bound on x for stage 0,',...
                    ' size [{{ dims_0.nbx_0 }}]\n ');
sfun_input_names = [sfun_input_names; 'lbx_0 [{{ dims_0.nbx_0 }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if dims_0.nbx_0 > 0 and simulink_opts.inputs.ubx_0 -%}  {#- ubx_0 #}
input_note = strcat(input_note, num2str(i_in), ') ubx_0 - upper bound on x for stage 0,',...
                    ' size [{{ dims_0.nbx_0 }}]\n ');
sfun_input_names = [sfun_input_names; 'ubx_0 [{{ dims_0.nbx_0 }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if np_total > 0 and simulink_opts.inputs.parameter_traj -%}  {#- parameter_traj #}
input_note = strcat(input_note, num2str(i_in), ') parameters - concatenated for all stages 0 to N,',...
                    ' size [{{ np_total }}]\n ');
sfun_input_names = [sfun_input_names; 'parameter_traj [{{ np_total }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if dims_0.np_global > 0 and simulink_opts.inputs.p_global -%}  {#- p_global #}
input_note = strcat(input_note, num2str(i_in), ') global parameters - first value indicates if update should be performed (0 means no update)\n');
input_note = strcat(input_note, '\tafterwards: new numerical values of p_global, size [1 + {{ dims_0.np_global }}]\n');
sfun_input_names = [sfun_input_names; 'p_global [1 + {{ dims_0.np_global }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if dims_0.ny_0 > 0 and simulink_opts.inputs.y_ref_0 %}
input_note = strcat(input_note, num2str(i_in), ') y_ref_0 - size [{{ dims_0.ny_0 }}]\n ');
sfun_input_names = [sfun_input_names; 'y_ref_0 [{{ dims_0.ny_0 }}]'];
i_in = i_in + 1;
{%- endif %}


{% if problem_class == "OCP" %}
{%- if dims.ny > 0 and solver_options.N_horizon > 1 and simulink_opts.inputs.y_ref %}
input_note = strcat(input_note, num2str(i_in), ') y_ref - concatenated for stages 1 to N-1,',...
                    ' size [{{ (solver_options.N_horizon-1) * dims.ny }}]\n ');
sfun_input_names = [sfun_input_names; 'y_ref [{{ (solver_options.N_horizon-1) * dims.ny }}]'];
i_in = i_in + 1;
{%- endif %}
{%- endif -%}

{%- if dims_e.ny_e > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.y_ref_e %}
input_note = strcat(input_note, num2str(i_in), ') y_ref_e - size [{{ dims_e.ny_e }}]\n ');
sfun_input_names = [sfun_input_names; 'y_ref_e [{{ dims_e.ny_e }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if nbx_total and simulink_opts.inputs.lbx -%}  {#- lbx #}
input_note = strcat(input_note, num2str(i_in), ') lbx values concatenated for stages 1 to N-1, size [{{ nbx_total }}]\n ');
sfun_input_names = [sfun_input_names; 'lbx [{{ nbx_total }}]'];
i_in = i_in + 1;
{%- endif %}
{%- if nbx_total and simulink_opts.inputs.ubx -%}  {#- ubx #}
input_note = strcat(input_note, num2str(i_in), ') ubx values concatenated for stages 1 to N-1, size [{{ nbx_total }}]\n ');
sfun_input_names = [sfun_input_names; 'ubx [{{ nbx_total }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if dims_e.nbx_e > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.lbx_e -%}  {#- lbx_e #}
input_note = strcat(input_note, num2str(i_in), ') lbx_e (lbx at shooting node N), size [{{ dims_e.nbx_e }}]\n ');
sfun_input_names = [sfun_input_names; 'lbx_e [{{ dims_e.nbx_e }}]'];
i_in = i_in + 1;
{%- endif %}
{%- if dims_e.nbx_e > 0 and solver_options.N_horizon > 0 and simulink_opts.inputs.ubx_e -%}  {#- ubx_e #}
input_note = strcat(input_note, num2str(i_in), ') ubx_e (ubx at shooting node N), size [{{ dims_e.nbx_e }}]\n ');
sfun_input_names = [sfun_input_names; 'ubx_e [{{ dims_e.nbx_e }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if nbu_total and simulink_opts.inputs.lbu -%}  {#- lbu #}
input_note = strcat(input_note, num2str(i_in), ') lbu for stages 0 to N-1, size [{{ nbu_total }}]\n ');
sfun_input_names = [sfun_input_names; 'lbu [{{ nbu_total }}]'];
i_in = i_in + 1;
{%- endif -%}
{%- if nbu_total and simulink_opts.inputs.ubu -%}  {#- ubu #}
input_note = strcat(input_note, num2str(i_in), ') ubu for stages 0 to N-1, size [{{ nbu_total }}]\n ');
sfun_input_names = [sfun_input_names; 'ubu [{{ nbu_total }}]'];
i_in = i_in + 1;
{%- endif -%}

{% if problem_class == "OCP" %}
{%- if dims.ng > 0 and simulink_opts.inputs.lg -%}  {#- lg #}
input_note = strcat(input_note, num2str(i_in), ') lg for stages 0 to N-1, size [{{ solver_options.N_horizon*dims.ng }}]\n ');
sfun_input_names = [sfun_input_names; 'lg [{{ solver_options.N_horizon*dims.ng }}]'];
i_in = i_in + 1;
{%- endif %}
{%- if dims.ng > 0 and simulink_opts.inputs.ug -%}  {#- ug #}
input_note = strcat(input_note, num2str(i_in), ') ug for stages 0 to N-1, size [{{ solver_options.N_horizon*dims.ng }}]\n ');
sfun_input_names = [sfun_input_names; 'ug [{{ solver_options.N_horizon*dims.ng }}]'];
i_in = i_in + 1;
{%- endif %}
{%- endif %}

{%- if nh_total > 0 and simulink_opts.inputs.lh -%}  {#- lh #}
input_note = strcat(input_note, num2str(i_in), ') lh for stages 1 to N-1, size [{{ nh_total }}]\n ');
sfun_input_names = [sfun_input_names; 'lh [{{ nh_total }}]'];
i_in = i_in + 1;
{%- endif %}
{%- if nh_total > 0 and simulink_opts.inputs.uh -%}  {#- uh #}
input_note = strcat(input_note, num2str(i_in), ') uh for stages 1 to N-1, size [{{ nh_total }}]\n ');
sfun_input_names = [sfun_input_names; 'uh [{{ nh_total }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if dims_0.nh_0 > 0 and simulink_opts.inputs.lh_0 -%}  {#- lh_0 #}
input_note = strcat(input_note, num2str(i_in), ') lh_0, size [{{ dims_0.nh_0 }}]\n ');
sfun_input_names = [sfun_input_names; 'lh_0 [{{ dims_0.nh_0 }}]'];
i_in = i_in + 1;
{%- endif %}
{%- if dims_0.nh_0 > 0 and simulink_opts.inputs.uh_0 -%}  {#- uh_0 #}
input_note = strcat(input_note, num2str(i_in), ') uh_0, size [{{ dims_0.nh_0 }}]\n ');
sfun_input_names = [sfun_input_names; 'uh_0 [{{ dims_0.nh_0 }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if dims_e.nh_e > 0 and simulink_opts.inputs.lh_e -%}  {#- lh_e #}
input_note = strcat(input_note, num2str(i_in), ') lh_e, size [{{ dims_e.nh_e }}]\n ');
sfun_input_names = [sfun_input_names; 'lh_e [{{ dims_e.nh_e }}]'];
i_in = i_in + 1;
{%- endif %}
{%- if dims_e.nh_e > 0 and simulink_opts.inputs.uh_e -%}  {#- uh_e #}
input_note = strcat(input_note, num2str(i_in), ') uh_e, size [{{ dims_e.nh_e }}]\n ');
sfun_input_names = [sfun_input_names; 'uh_e [{{ dims_e.nh_e }}]'];
i_in = i_in + 1;
{%- endif %}


{% if problem_class == "OCP" %}
{%- if dims_0.ny_0 > 0 and simulink_opts.inputs.cost_W_0 %}  {#- cost_W_0 #}
input_note = strcat(input_note, num2str(i_in), ') cost_W_0 in column-major format, size [{{ dims_0.ny_0 * dims_0.ny_0 }}]\n ');
sfun_input_names = [sfun_input_names; 'cost_W_0 [{{ dims_0.ny_0 * dims_0.ny_0 }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if dims.ny > 0 and simulink_opts.inputs.cost_W %}  {#- cost_W #}
input_note = strcat(input_note, num2str(i_in), ') cost_W in column-major format, that is set for all intermediate stages: 1 to N-1, size [{{ dims.ny * dims.ny }}]\n ');
sfun_input_names = [sfun_input_names; 'cost_W [{{ dims.ny * dims.ny }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if dims_e.ny_e > 0 and simulink_opts.inputs.cost_W_e %}  {#- cost_W_e #}
input_note = strcat(input_note, num2str(i_in), ') cost_W_e in column-major format, size [{{ dims_e.ny_e * dims_e.ny_e }}]\n ');
sfun_input_names = [sfun_input_names; 'cost_W_e [{{ dims_e.ny_e * dims_e.ny_e }}]'];
i_in = i_in + 1;
{%- endif %}
{%- endif %}


{%- if ns_total > 0 and simulink_opts.inputs.cost_zl %}  {#- cost_zl #}
input_note = strcat(input_note, num2str(i_in), ') cost_zl for all nodes 0 to N, size [{{ ns_total }}]\n ');
sfun_input_names = [sfun_input_names; 'cost_zl [{{ ns_total }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if ns_total > 0 and simulink_opts.inputs.cost_zu %}  {#- cost_zu #}
input_note = strcat(input_note, num2str(i_in), ') cost_zu for all nodes 0 to N, size [{{ ns_total }}]\n ');
sfun_input_names = [sfun_input_names; 'cost_zu [{{ ns_total }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if ns_total > 0 and simulink_opts.inputs.cost_Zl %}  {#- cost_Zl #}
input_note = strcat(input_note, num2str(i_in), ') cost_Zl for all nodes 0 to N, size [{{ ns_total }}]\n ');
sfun_input_names = [sfun_input_names; 'cost_Zl [{{ ns_total }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if ns_total > 0 and simulink_opts.inputs.cost_Zu %}  {#- cost_Zu #}
input_note = strcat(input_note, num2str(i_in), ') cost_Zu for all nodes 0 to N, size [{{ ns_total }}]\n ');
sfun_input_names = [sfun_input_names; 'cost_Zu [{{ ns_total }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if simulink_opts.inputs.reset_solver %}  {#- reset_solver #}
input_note = strcat(input_note, num2str(i_in), ') reset_solver - determines if iterate is set to all zeros before other initializations (x_init, u_init, pi_init) are set and before solver is called, size [1]\n ');
sfun_input_names = [sfun_input_names; 'reset_solver [1]'];
i_in = i_in + 1;
{%- endif %}

{%- if simulink_opts.inputs.ignore_inits %}  {#- ignore_inits #}
input_note = strcat(input_note, num2str(i_in), ') ignore_inits - determines if initialization (x_init, u_init, pi_init, slacks_init) are set (ignore_inits == 0) or ignored (otherwise), ignoring corresponds to internal warm start, size [1]\n ');
sfun_input_names = [sfun_input_names; 'ignore_inits [1]'];
i_in = i_in + 1;
{%- endif %}

{%- if simulink_opts.inputs.x_init %}  {#- x_init #}
input_note = strcat(input_note, num2str(i_in), ') x_init - initialization of x for all stages, size [{{ nx_total }}]\n ');
sfun_input_names = [sfun_input_names; 'x_init [{{ nx_total }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if simulink_opts.inputs.u_init %}  {#- u_init #}
input_note = strcat(input_note, num2str(i_in), ') u_init - initialization of u for stages 0 to N-1, size [{{ nu_total }}]\n ');
sfun_input_names = [sfun_input_names; 'u_init [{{ nu_total }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if simulink_opts.inputs.pi_init %}  {#- pi_init #}
input_note = strcat(input_note, num2str(i_in), ') pi_init - initialization of pi for stages 0 to N-1, size [{{ npi_total }}]\n ');
sfun_input_names = [sfun_input_names; 'pi_init [{{ npi_total }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if simulink_opts.inputs.slacks_init %}  {#- slacks_init #}
input_note = strcat(input_note, num2str(i_in), ') slacks_init - initialization of slack values for all stages (0 to N), size [{{ two_ns_total }}]');
sfun_input_names = [sfun_input_names; 'slacks_init [{{ two_ns_total }}]'];
i_in = i_in + 1;
{%- endif %}

{%- if simulink_opts.inputs.rti_phase %}  {#- rti_phase #}
input_note = strcat(input_note, num2str(i_in), ') rti_phase, size [1]\n ');
sfun_input_names = [sfun_input_names; 'rti_phase [1]'];
i_in = i_in + 1;
{%- endif %}

{%- if simulink_opts.inputs.levenberg_marquardt %}  {#- levenberg_marquardt #}
input_note = strcat(input_note, num2str(i_in), ') levenberg_marquardt, size [1]\n ');
sfun_input_names = [sfun_input_names; 'levenberg_marquardt [1]'];
i_in = i_in + 1;
{%- endif %}

{%- if simulink_opts.customizable_inputs %}
{#- customizable inputs #}
{%- for input_name, input_spec in simulink_opts.customizable_inputs -%}
  {%- if input_name is starting_with("sparse_parameter_") %}
    {% set param_length = input_spec.parameter_indices | length %}
    {% set port_name = input_name | replace(from="sparse_parameter_", to="") %}
    {% set stage_idx_0 = input_spec.stage_idx_0 %}
    {% set stage_idx_e = input_spec.stage_idx_e %}
    % {{ port_name }}
input_note = strcat(input_note, num2str(i_in), ') {{ input_name }}, to update parameter values at nodes {{ stage_idx_0 }} to {{ stage_idx_e }}, \n\twith {{ param_length }} parameter indices specified [{{ input_spec.parameter_indices | first }}, ..., {{ input_spec.parameter_indices | last }}]\n ');
input_note = strcat(input_note, '\tfirst value indicates if update should be performed (0 means no update)\n');
input_note = strcat(input_note, '\tafterwards: new numerical values of parameters to be updated at stages, size [1 + {{ (stage_idx_e - stage_idx_0 + 1) }} * {{ param_length }} = {{ 1 + (stage_idx_e - stage_idx_0 + 1) * param_length }}]\n');
sfun_input_names = [sfun_input_names; 'sparse_param_{{ port_name }}'];
i_in = i_in + 1;
  {%- endif -%}
{%- endfor %}
{%- endif -%}

fprintf(input_note)

disp(' ')

output_note = 'Outputs are:\n';
i_out = 0;

global sfun_output_names
sfun_output_names = {};

{%- if dims_0.nu > 0 and simulink_opts.outputs.u0 == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') u0, control input at node 0, size [{{ dims_0.nu }}]\n ');
sfun_output_names = [sfun_output_names; 'u0 [{{ dims_0.nu }}]'];
{%- endif %}

{%- if simulink_opts.outputs.utraj == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') utraj, control input concatenated for nodes 0 to N-1, size [{{ nu_total }}]\n ');
sfun_output_names = [sfun_output_names; 'utraj [{{ nu_total }}]'];
{%- endif %}

{%- if simulink_opts.outputs.xtraj == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') xtraj, state concatenated for nodes 0 to N, size [{{ nx_total }}]\n ');
sfun_output_names = [sfun_output_names; 'xtraj [{{ nx_total }}]'];
{%- endif %}

{%- if simulink_opts.outputs.ztraj == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') ztraj, algebraic states concatenated for nodes 0 to N-1, size [{{ nz_total }}]\n ');
sfun_output_names = [sfun_output_names; 'ztraj [{{ nz_total }}]'];
{%- endif %}

{%- if simulink_opts.outputs.pi_all == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') pi_all, equality Lagrange multipliers concatenated for nodes 0 to N-1, size [{{ npi_total }}]\n ');
sfun_output_names = [sfun_output_names; 'pi_all [{{ npi_total }}]'];
{%- endif %}

{%- if simulink_opts.outputs.slack_values == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') slack values concatenated in order [sl_0, su_0, ..., sl_N, su_N] \n ');
sfun_output_names = [sfun_output_names; 'slack_values [{{ two_ns_total }}]'];
{%- endif %}

{%- if simulink_opts.outputs.solver_status == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') acados solver status (0 = SUCCESS)\n ');
sfun_output_names = [sfun_output_names; 'solver_status'];
{%- endif %}

{%- if simulink_opts.outputs.cost_value == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') cost function value\n ');
sfun_output_names = [sfun_output_names; 'cost_value'];
{%- endif %}


{%- if simulink_opts.outputs.KKT_residual == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') KKT residual\n ');
sfun_output_names = [sfun_output_names; 'KKT_residual'];
{%- endif %}

{%- if simulink_opts.outputs.KKT_residuals == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') KKT residuals, size [4] (stat, eq, ineq, comp)\n ');
sfun_output_names = [sfun_output_names; 'KKT_residuals [4]'];
{%- endif %}

{%- if solver_options.N_horizon > 0 and simulink_opts.outputs.x1 == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') x1, state at node 1\n ');
sfun_output_names = [sfun_output_names; 'x1 [{{ dims_0.nx_next }}]'];
{%- endif %}

{%- if simulink_opts.outputs.CPU_time == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') CPU time\n ');
sfun_output_names = [sfun_output_names; 'CPU_time'];
{%- endif %}

{%- if simulink_opts.outputs.CPU_time_sim == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') CPU time integrator\n ');
sfun_output_names = [sfun_output_names; 'CPU_time_sim'];
{%- endif %}

{%- if simulink_opts.outputs.CPU_time_qp == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') CPU time QP solution\n ');
sfun_output_names = [sfun_output_names; 'CPU_time_qp'];
{%- endif %}

{%- if simulink_opts.outputs.CPU_time_lin == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') CPU time linearization (including integrator)\n ');
sfun_output_names = [sfun_output_names; 'CPU_time_lin'];
{%- endif %}

{%- if simulink_opts.outputs.sqp_iter == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') SQP iterations\n ');
sfun_output_names = [sfun_output_names; 'sqp_iter'];
{%- endif %}

{%- if simulink_opts.outputs.parameter_traj == 1 %}
i_out = i_out + 1;
output_note = strcat(output_note, num2str(i_out), ') parameter trajectory\n ');
sfun_output_names = [sfun_output_names; 'parameter_traj [{{ np_total }}]'];
{%- endif %}

fprintf(output_note)

{%- if simulink_opts.generate_simulink_block == 1 %}
modelName = '{{ name }}_ocp_solver_simulink_block';
new_system(modelName);
open_system(modelName);

blockPath = [modelName '/{{ name }}_ocp_solver'];
add_block('simulink/User-Defined Functions/S-Function', blockPath);
set_param(blockPath, 'FunctionName', 'acados_solver_sfunction_{{ name }}');

Simulink.Mask.create(blockPath);
{%- if simulink_opts.show_port_info == 1 %}
mask_str = sprintf([ ...
    'global sfun_input_names sfun_output_names\n' ...
    'for i = 1:length(sfun_input_names)\n' ...
    '    port_label(''input'', i, sfun_input_names{i})\n' ...
    'end\n' ...
    'for i = 1:length(sfun_output_names)\n' ...
    '    port_label(''output'', i, sfun_output_names{i})\n' ...
    'end\n' ...
    'disp("acados OCP")' ...
]);
{%- else %}
mask_str = sprintf('disp("acados OCP")');
{%- endif %}
mask = Simulink.Mask.get(blockPath);
mask.Display = mask_str;

save_system(modelName);
close_system(modelName);
disp([newline, 'Created the OCP solver Simulink block in: ', modelName])
{%- endif %}

% The mask drawing command is:
% ---
% global sfun_input_names sfun_output_names
% for i = 1:length(sfun_input_names)
%     port_label('input', i, sfun_input_names{i})
% end
% for i = 1:length(sfun_output_names)
%     port_label('output', i, sfun_output_names{i})
% end
% ---
% It can be used by copying it in sfunction/Mask/Edit mask/Icon drawing commands
%   (you can access it with ctrl+M on the s-function)
