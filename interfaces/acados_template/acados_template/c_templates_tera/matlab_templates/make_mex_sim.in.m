%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

function make_mex_sim_{{ name }}()

    opts.output_dir = pwd;

    % get acados folder
    acados_folder = getenv('ACADOS_INSTALL_DIR');

    % set paths
    acados_include = ['-I' fullfile(acados_folder, 'include')];
    template_lib_include = ['-l' 'acados_sim_solver_{{ name }}'];
    template_lib_path = ['-L' fullfile(pwd)];

    acados_link_str = ['-L' '{{ code_gen_options.acados_lib_path }}'];
    external_include = ['-I', fullfile(acados_folder,'include')];
    blasfeo_include = ['-I', fullfile(acados_folder,'include', 'blasfeo', 'include')];
    hpipm_include = ['-I', fullfile(acados_folder,'include', 'hpipm', 'include')];

    % load linking information of compiled acados
    link_libs_core_filename = fullfile(acados_folder, 'lib', 'link_libs.json');
    addpath(fullfile(acados_folder, 'external', 'jsonlab'));
    link_libs = loadjson(link_libs_core_filename);

    % add necessary link instructs
    acados_def_extra = '';
    acados_lib_extra = {};
    lib_names = fieldnames(link_libs);
    for idx = 1 : numel(lib_names)
        lib_name = lib_names{idx};
        link_arg = link_libs.(lib_name);
        if ~isempty(link_arg)
            def_arg = sprintf('-DACADOS_WITH_%s', upper(lib_name));
            acados_def_extra = [acados_def_extra, ' ', def_arg];
            if ~strcmp(lib_name, 'openmp')
                acados_lib_extra = [acados_lib_extra, link_arg];
            end
        end
    end

    mex_include = ['-I', fullfile(acados_folder, 'interfaces', 'acados_matlab_octave')];

    mex_names = { ...
        'acados_sim_create_{{ name }}', ...
        'acados_sim_free_{{ name }}', ...
        'acados_sim_set_{{ name }}' ...
    };

    mex_files = cell(length(mex_names), 1);
    for k=1:length(mex_names)
        mex_files{k} = fullfile([mex_names{k}, '.c']);
    end

    %% octave C flags
    if is_octave()
        if ~exist(fullfile(opts.output_dir, 'cflags_octave.txt'), 'file')
            diary(fullfile(opts.output_dir, 'cflags_octave.txt'));
            diary on
            mkoctfile -p CFLAGS
            diary off
            input_file = fopen(fullfile(opts.output_dir, 'cflags_octave.txt'), 'r');
            cflags_tmp = fscanf(input_file, '%[^\n]s');
            fclose(input_file);
            if ~ismac()
                cflags_tmp = [cflags_tmp, ' -std=c99 -fopenmp'];
            else
                cflags_tmp = [cflags_tmp, ' -std=c99'];
            end
            input_file = fopen(fullfile(opts.output_dir, 'cflags_octave.txt'), 'w');
            fprintf(input_file, '%s', cflags_tmp);
            fclose(input_file);
        end
        % read cflags from file
        input_file = fopen(fullfile(opts.output_dir, 'cflags_octave.txt'), 'r');
        cflags_tmp = fscanf(input_file, '%[^\n]s');
        fclose(input_file);
        setenv('CFLAGS', cflags_tmp);
    end

    %% flags
    FLAGS = 'CFLAGS=$CFLAGS -std=c99';
    LDFLAGS = 'LDFLAGS=$LDFLAGS';
    COMPFLAGS = 'COMPFLAGS=$COMPFLAGS';
    COMPDEFINES = 'COMPDEFINES=$COMPDEFINES';
    if ~ismac() && ~isempty(link_libs.openmp)
        LDFLAGS = [LDFLAGS, ' ', link_libs.openmp];
        COMPFLAGS = [COMPFLAGS, ' ', link_libs.openmp]; % seems unnecessary
    end
    if ~is_octave()
        FLAGS = [FLAGS, acados_def_extra];
        COMPDEFINES = [COMPDEFINES, acados_def_extra];
    end

    %% compile mex
    for ii=1:length(mex_files)
        disp(['compiling ', mex_files{ii}])
        if is_octave()
    %        mkoctfile -p CFLAGS
            mex(acados_include, template_lib_include, external_include, blasfeo_include, hpipm_include,...
                template_lib_path, mex_include, acados_link_str, '-lacados', '-lhpipm', '-lblasfeo',...
                acados_lib_extra{:}, mex_files{ii})
        else
            mex(FLAGS, LDFLAGS, COMPDEFINES, COMPFLAGS, acados_include, template_lib_include, external_include, blasfeo_include, hpipm_include,...
                template_lib_path, mex_include, acados_link_str, '-lacados', '-lhpipm', '-lblasfeo',...
                acados_lib_extra{:}, mex_files{ii})
        end
    end
end
