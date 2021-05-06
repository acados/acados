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

function ocp_compile_interface(opts)

% get acados folder
acados_folder = getenv('ACADOS_INSTALL_DIR');
mex_flags = getenv('ACADOS_MEX_FLAGS');

% set paths
acados_mex_folder = fullfile(acados_folder, 'interfaces', 'acados_matlab_octave');
acados_include = ['-I', acados_folder];
acados_interfaces_include = ['-I', fullfile(acados_folder, 'interfaces')];
external_include = ['-I', fullfile(acados_folder, 'external')];
blasfeo_include = ['-I', fullfile(acados_folder, 'external', 'blasfeo', 'include')];
hpipm_include = ['-I', fullfile(acados_folder, 'external', 'hpipm', 'include')];
acados_lib_path = ['-L', fullfile(acados_folder, 'lib')];

mex_names = { ...
    'ocp_create', ...
    'ocp_destroy', ...
    'ocp_create_ext_fun', ...
    'ocp_destroy_ext_fun', ...
    'ocp_solve', ...
    'ocp_get_cost', ...
    'ocp_precompute', ...
    'ocp_set', ...
    'ocp_get' ...
    'ocp_eval_param_sens', ...
};
mex_files = cell(length(mex_names), 1);
for k=1:length(mex_names)
    mex_files{k} = fullfile(acados_mex_folder, [mex_names{k}, '.c']);
end

%% check linking information of compiled acados
% copy link_libs.json to build to check for consistency later
link_libs_core_filename = fullfile(acados_folder, 'lib', 'link_libs.json');
link_libs_interface_filename = fullfile(opts.output_dir, 'link_libs.json');
copyfile(link_libs_core_filename, link_libs_interface_filename);
addpath(fullfile(acados_folder, 'external', 'jsonlab'));
libs = loadjson(link_libs_core_filename);

%% compile mex
if is_octave()
    if ~exist(fullfile(opts.output_dir, 'cflags_octave.txt'), 'file')
        diary(fullfile(opts.output_dir, 'cflags_octave.txt'));
        diary on
        mkoctfile -p CFLAGS
        diary off
        input_file = fopen(fullfile(opts.output_dir, 'cflags_octave.txt'), 'r');
        cflags_tmp = fscanf(input_file, '%[^\n]s');
        fclose(input_file);
        cflags_tmp = [cflags_tmp, ' -std=c99'];
        input_file = fopen(fullfile(opts.output_dir, 'cflags_octave.txt'), 'w');
        fprintf(input_file, '%s', cflags_tmp);
        fclose(input_file);
    end
    % read cflags from file
    input_file = fopen(fullfile(opts.output_dir, 'cflags_octave.txt'), 'r');
    cflags_tmp = fscanf(input_file, '%[^\n]s');
    fclose(input_file);

    % add flags
    defines_tmp = '';
    if ~isempty(libs.qpoases)
        defines_tmp = [defines_tmp, ' -DACADOS_WITH_QPOASES'];
    end
    if ~isempty(libs.osqp)
        defines_tmp = [defines_tmp, ' -DACADOS_WITH_OSQP'];
    end
    if ~isempty(libs.hpmpc)
        defines_tmp = [defines_tmp, ' -DACADOS_WITH_HPMPC'];
    end
    if ~isempty(libs.qpdunes)
        defines_tmp = [defines_tmp, ' -DACADOS_WITH_QPDUNES'];
    end
    if ~isempty(libs.ooqp)
        defines_tmp = [defines_tmp, ' -DACADOS_WITH_OOQP'];
    end
    setenv('CFLAGS', [cflags_tmp, defines_tmp]);
    setenv('COMPDEFINES', defines_tmp);

    if ~ismac() && ~isempty(libs.openmp)
        setenv('LDFLAGS', libs.openmp);
        setenv('COMPFLAGS', libs.openmp);
    end

end


FLAGS = 'CFLAGS=$CFLAGS -std=c99';
LDFLAGS = 'LDFLAGS=$LDFLAGS';
COMPFLAGS = 'COMPFLAGS=$COMPFLAGS';
COMPDEFINES = 'COMPDEFINES=$COMPDEFINES';
if ~ismac() && ~isempty(libs.openmp)
    LDFLAGS = [LDFLAGS, ' ', libs.openmp];
    COMPFLAGS = [COMPFLAGS, ' ', libs.openmp]; % seems unnecessary
end

if ~is_octave()
    defines_tmp = '';
    if ~isempty(libs.qpoases)
        defines_tmp = [defines_tmp, ' -DACADOS_WITH_QPOASES'];
    end
    if ~isempty(libs.osqp)
        defines_tmp = [defines_tmp, ' -DACADOS_WITH_OSQP'];
    end
    if ~isempty(libs.hpmpc)
        defines_tmp = [defines_tmp, ' -DACADOS_WITH_HPMPC'];
    end
    if ~isempty(libs.qpdunes)
        defines_tmp = [defines_tmp, ' -DACADOS_WITH_QPDUNES'];
    end
    if ~isempty(libs.ooqp)
        defines_tmp = [defines_tmp, ' -DACADOS_WITH_OOQP'];
    end
    FLAGS = [FLAGS, defines_tmp];
    COMPDEFINES = [COMPDEFINES, defines_tmp];
end

for ii=1:length(mex_files)
    disp(['compiling ', mex_files{ii}])
    if is_octave()
        fn = fieldnames(libs);
        linker_flags = {'-lacados', '-lhpipm', '-lblasfeo'};
        for k = 1:numel(fn)
            if ~isempty(libs.(fn{k}))
                linker_flags{end+1} = libs.(fn{k});
            end
        end
        % NOTE: multiple linker flags in 1 argument do not work in Matlab
        mex(acados_include, acados_interfaces_include, external_include, blasfeo_include, hpipm_include,...
            acados_lib_path, linker_flags{:}, mex_files{ii})
    else
        % gcc uses FLAGS, LDFLAGS
        % MSVC uses COMPFLAGS, COMPDEFINES
        % NOTE: empty linker flags do not work in Octave
        mex(mex_flags, FLAGS, LDFLAGS, COMPDEFINES, COMPFLAGS, acados_include, acados_interfaces_include, external_include, blasfeo_include, hpipm_include,...
            acados_lib_path, '-lacados', '-lhpipm', '-lblasfeo', libs.qpoases,...
            libs.qpdunes, libs.osqp, libs.hpmpc, libs.ooqp, mex_files{ii}, '-outdir', opts.output_dir)
    end
end

if is_octave()
    octave_version = OCTAVE_VERSION();
    if octave_version < 5
        movefile('*.o', opts.output_dir);
    end

    %system(['mv -f *.mexa64 ', opts.output_dir])
    for k=1:length(mex_names)
        clear(mex_names{k})
    %    movefile([mex_names{k}, '.', mexext], opts.output_dir);
        [status, message] = copyfile([mex_names{k}, '.', mexext], opts.output_dir);
        delete([mex_names{k}, '.', mexext]);
    end
end


