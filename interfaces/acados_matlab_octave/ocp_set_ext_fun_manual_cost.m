%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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

function C_ocp_ext_fun = ocp_set_ext_fun_manual_cost(C_ocp, C_ocp_ext_fun, model_struct, opts_struct, ...
    ext_fun_type, cost_file_name, cost_function_name, cost_e_file_name, cost_e_function_name)

model_name = model_struct.name;
N = opts_struct.param_scheme_N;

% get acados folder
acados_folder = getenv('ACADOS_INSTALL_DIR');
mex_flags = getenv('ACADOS_MEX_FLAGS');

% set paths
acados_mex_folder = fullfile(acados_folder, 'interfaces', 'acados_matlab_octave');
acados_include = ['-I' acados_folder];
acados_interfaces_include = ['-I' fullfile(acados_folder, 'interfaces')];
external_include = ['-I' fullfile(acados_folder, 'external')];
blasfeo_include = ['-I' fullfile(acados_folder, 'external' , 'blasfeo', 'include')];
hpipm_include = ['-I' fullfile(acados_folder, 'external' , 'hpipm', 'include')];
acados_lib_path = ['-L' fullfile(acados_folder, 'lib')];
acados_matlab_octave_lib_path = ['-L' fullfile(acados_folder, 'interfaces', 'acados_matlab_octave')];
model_lib_path = ['-L', opts_struct.output_dir];

%% select files to compile
%mex_files = {};
% ext_fun_type = 'casadi' or generic
mex_files = {fullfile(acados_mex_folder, ['ocp_set_ext_fun_', ext_fun_type, '.c'])};
setter = {};
set_fields = {};
mex_fields = {};
fun_names = {};
mex_names = {};
source_file_names = {};
phase = {};
phase_start = {};
phase_end = {};

% external cost
if (strcmp(model_struct.cost_type, 'ext_cost'))

    setter = {setter{:} ...        
        'ocp_nlp_cost_model_set' ...        
        };
    set_fields = {set_fields{:} ...        
        'ext_cost_fun_jac_hess' ...        
        };
    mex_fields = {mex_fields{:} ...        
        'cost_ext_cost_fun_jac_hess' ...        
        };
    fun_names = {fun_names{:} ...        
        cost_function_name ... %[model_name, '_cost_ext_cost_fun_jac_hess'] ...        
        };
    mex_names = {mex_names{:} ...       
        [model_name, '_ocp_set_ext_fun_cost_0_fun_jac_hess'] ...        
        };
    source_file_names = {source_file_names{:} ...        
        cost_file_name ...
        };    
    phase = {phase{:}, 0};
    phase_start = {phase_start{:}, 0};
    phase_end = {phase_end{:}, N-1};

end
if (strcmp(model_struct.cost_type_e, 'ext_cost'))

    setter = {setter{:} ...        
        'ocp_nlp_cost_model_set' ...
        %'ocp_nlp_cost_model_set' ...
        };
    set_fields = {set_fields{:} ...        
        'ext_cost_fun_jac_hess' ...
        %'ext_cost_fun' ...
        };
    mex_fields = {mex_fields{:} ...        
        'cost_ext_cost_fun_jac_hess' ...
        %'cost_ext_cost_fun' ...
        };
    fun_names = {fun_names{:} ...        
        cost_e_function_name ... %[model_name, '_cost_ext_cost_e_fun_jac_hess'] ...
        %[model_name, '_cost_ext_cost_e_fun'] ...
        };
    mex_names = {mex_names{:} ...        
        [model_name, '_ocp_set_ext_fun_cost_1_fun_jac_hess'] ...
        %[model_name, '_ocp_set_ext_fun_cost_1_ext_cost_fun'] ...
        };
    source_file_names = {source_file_names{:} ...        
        cost_e_file_name ...
        };     
    phase = {phase{:}, 1};
    phase_start = {phase_start{:}, N};
    phase_end = {phase_end{:}, N};

end

% compile mex files
%mex_files
%setter
%set_fields
%mex_fields
%fun_names
%mex_names
%phase
%phase_start
%phase_end
if (strcmp(opts_struct.compile_interface, 'true') || strcmp(opts_struct.codgen_model, 'true'))

    if is_octave()
        if ~exist(fullfile(opts_struct.output_dir, 'cflags_octave.txt'), 'file')
            diary(fullfile(opts_struct.output_dir, 'cflags_octave.txt'))
            diary on
            mkoctfile -p CFLAGS
            diary off
            input_file = fopen(fullfile(opts_struct.output_dir, 'cflags_octave.txt'), 'r');
            cflags_tmp = fscanf(input_file, '%[^\n]s');
            fclose(input_file);
            if ~ismac()
                cflags_tmp = [cflags_tmp, ' -std=c99 -fopenmp'];
            else
                cflags_tmp = [cflags_tmp, ' -std=c99'];
            end
            input_file = fopen(fullfile(opts_struct.output_dir, 'cflags_octave.txt'), 'w');
            fprintf(input_file, '%s', cflags_tmp);
            fclose(input_file);
        end
    end

    %% get pointers for external functions in model
    for ii=1:length(mex_names)

        disp(['compiling ', mex_names{ii}])
        if is_octave()
    %        mkoctfile -p CFLAGS
            input_file = fopen(fullfile(opts_struct.output_dir, 'cflags_octave.txt'), 'r');
            cflags_tmp = fscanf(input_file, '%[^\n]s');
            fclose(input_file);
            cflags_tmp = [cflags_tmp, ' -DSETTER=', setter{ii}];
            cflags_tmp = [cflags_tmp, ' -DSET_FIELD=', set_fields{ii}];
            cflags_tmp = [cflags_tmp, ' -DMEX_FIELD=', mex_fields{ii}];
            cflags_tmp = [cflags_tmp, ' -DFUN_NAME=', fun_names{ii}];
            cflags_tmp = [cflags_tmp, ' -DPHASE=', num2str(phase{ii})];
            cflags_tmp = [cflags_tmp, ' -DN0=', num2str(phase_start{ii})];
            cflags_tmp = [cflags_tmp, ' -DN1=', num2str(phase_end{ii})];
            setenv('CFLAGS', cflags_tmp);
            mex(acados_include, acados_interfaces_include, external_include, blasfeo_include,...
                hpipm_include, acados_lib_path, acados_matlab_octave_lib_path, model_lib_path, '-lacados',...
                 '-lhpipm', '-lblasfeo', source_file_names{ii}, mex_files{1});
        else
            if ~ismac()
                FLAGS = 'CFLAGS=$CFLAGS -std=c99 -fopenmp';
            else
                FLAGS = 'CFLAGS=$CFLAGS -std=c99 -O2 -fPIC';
            end
            mex(mex_flags, FLAGS, ['-DSETTER=', setter{ii}],...
                ['-DSET_FIELD=', set_fields{ii}], ['-DMEX_FIELD=', mex_fields{ii}],...
                ['-DFUN_NAME=', fun_names{ii}], ['-DPHASE=', num2str(phase{ii})],...
                ['-DN0=', num2str(phase_start{ii})], ['-DN1=', num2str(phase_end{ii})],...
                acados_include, acados_interfaces_include, external_include, blasfeo_include,...
                hpipm_include, acados_lib_path, acados_matlab_octave_lib_path, model_lib_path,...
                '-lacados', '-lhpipm', '-lblasfeo', source_file_names{ii}, mex_files{1});
            disable_last_warning();
        end
        
%        clear(mex_names{ii})
        [filepath,name,ext] = fileparts(source_file_names{ii});
        out_name = [name '.', mexext];
        %movefile(['ocp_set_ext_fun_', ext_fun_type, '.', mexext], fullfile(opts_struct.output_dir, [mex_names{ii}, '.', mexext]));
        movefile(out_name, fullfile(opts_struct.output_dir, ['manual_' mex_names{ii}, '.', mexext]));
        
    end
    
    if is_octave()
        octave_version = OCTAVE_VERSION();
        if octave_version < 5
            movefile('*.o', opts_struct.output_dir);
        end
    end

end
%C_ocp_ext_fun

% codegen the file to call mex files
%fileID = fopen('build/ocp_set_ext_fun_tmp.m', 'w');
%fprintf(fileID, 'function C_ocp_ext_fun = ocp_set_ext_fun_tmp(C_ocp, C_ocp_ext_fun, model_struct, opts_struct)\n');
for ii=1:length(mex_names)
%     disp(['evaluating: ' mex_names{ii}, '(C_ocp, C_ocp_ext_fun, model_struct, opts_struct)'])
    C_ocp_ext_fun = eval(['manual_' mex_names{ii} '(C_ocp, C_ocp_ext_fun, model_struct, opts_struct)']);
%    disp(['eval ', mex_names{ii}, ' done']);
%    fprintf(fileID, [mex_names{ii}, '(C_ocp, C_ocp_ext_fun, model_struct, opts_struct);\n']);
end
%fclose(fileID);


end
