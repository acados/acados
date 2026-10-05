%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%


{%- if solver_options.hessian_approx %}
    {%- set hessian_approx = solver_options.hessian_approx %}
{%- elif solver_options.sens_hess %}
    {%- set hessian_approx = "EXACT" %}
{%- else %}
    {%- set hessian_approx = "GAUSS_NEWTON" %}
{%- endif %}

SOURCES = { ...
        'acados_sim_solver_sfunction_{{ name }}.c', ...
        'acados_sim_solver_{{ name }}.c', ...
{%- for filename in external_function_files_model %}
        '{{ filename }}', ...
{%- endfor %}
      };

INC_PATH = '{{ code_gen_options.acados_include_path }}';

INCS = {['-I', fullfile(INC_PATH, 'blasfeo', 'include')], ...
    ['-I', fullfile(INC_PATH, 'hpipm', 'include')], ...
    ['-I', fullfile(INC_PATH, 'acados')], ...
    ['-I', fullfile(INC_PATH)]};

CFLAGS = 'CFLAGS=$CFLAGS';
LDFLAGS = 'LDFLAGS=$LDFLAGS';
COMPFLAGS = 'COMPFLAGS=$COMPFLAGS';
COMPDEFINES = 'COMPDEFINES=$COMPDEFINES';

LIB_PATH = ['-L', fullfile('{{ code_gen_options.acados_lib_path }}')];

LIBS = {'-lacados', '-lhpipm', '-lblasfeo'};

COMPFLAGS = [COMPFLAGS ' {{ code_gen_options.ext_fun_compile_flags }}'];
CFLAGS = [CFLAGS ' {{ code_gen_options.ext_fun_compile_flags }}'];


try
    % mex('-v', '-O', CFLAGS, LDFLAGS, COMPFLAGS, COMPDEFINES, INCS{:}, ...
    mex('-O', CFLAGS, LDFLAGS, COMPFLAGS, COMPDEFINES, INCS{:}, ...
        LIB_PATH, LIBS{:}, SOURCES{:}, ...
        '-output', 'acados_sim_solver_sfunction_{{ name }}');

catch exception
    disp('make_sfun_sim failed with the following exception:')
    disp(exception);
    disp(exception.message);
    disp('Try adding -v to the mex command above to get more information.')
    keyboard
end


fprintf( [ '\n\nSuccessfully created sfunction:\nacados_sim_solver_sfunction_{{ name }}', '.', ...
    eval('mexext')] );


global sfun_sim_input_names
sfun_sim_input_names = {};

%% print note on usage of s-function
fprintf('\n\nNote: Usage of Sfunction is as follows:\n')
input_note = 'Inputs are:\n1) x0, initial state, size [{{ dims.nx }}]\n ';
i_in = 2;
sfun_sim_input_names = [sfun_sim_input_names; 'x0 [{{ dims.nx }}]'];

{%- if dims.nu > 0 %}
input_note = strcat(input_note, num2str(i_in), ') u, size [{{ dims.nu }}]\n ');
i_in = i_in + 1;
sfun_sim_input_names = [sfun_sim_input_names; 'u [{{ dims.nu }}]'];
{%- endif %}

{%- if dims.np > 0 %}
input_note = strcat(input_note, num2str(i_in), ') parameters, size [{{ dims.np }}]\n ');
i_in = i_in + 1;
sfun_sim_input_names = [sfun_sim_input_names; 'p [{{ dims.np }}]'];
{%- endif %}

fprintf(input_note)

disp(' ')

global sfun_sim_output_names
sfun_sim_output_names = {};

output_note = strcat('Outputs are:\n', ...
                '1) x1 - simulated state, size [{{ dims.nx }}]\n');
sfun_sim_output_names = [sfun_sim_output_names; 'x1 [{{ dims.nx }}]'];

fprintf(output_note)


% create the Simulink block for the integrator
modelName = '{{ name }}_sim_solver_simulink_block';
new_system(modelName);
open_system(modelName);

blockPath = [modelName '/{{ name }}_sim_solver'];
add_block('simulink/User-Defined Functions/S-Function', blockPath);
set_param(blockPath, 'FunctionName', 'acados_sim_solver_sfunction_{{ name }}');

Simulink.Mask.create(blockPath);


display_name = '{{ name }} acados sim';
input_labels = '';
for i = 1:length(sfun_sim_input_names)
	input_labels = [input_labels, sprintf('port_label(''input'', %d, ''%s'')\n', i, sfun_sim_input_names{i})];
end
output_labels = '';
for i = 1:length(sfun_sim_output_names)
	output_labels = [output_labels, sprintf('port_label(''output'', %d, ''%s'')\n', i, sfun_sim_output_names{i})];
end
mask_str = [input_labels, output_labels, sprintf('disp(''%s'')', display_name)];

mask = Simulink.Mask.get(blockPath);
mask.Display = mask_str;

save_system(modelName);
close_system(modelName);
disp([newline, 'Created the sim solver Simulink block in: ', modelName])
