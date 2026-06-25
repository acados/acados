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

classdef AcadosCodeGenOptions < handle

    properties
        acados_include_path
        acados_link_libs
        shared_lib_ext
        os
        cython_include_dirs
        acados_lib_path
        json_file
        code_export_directory
        acados_version
        casadi_code_gen_options

        ext_fun_compile_flags
        ext_fun_expand_constr
        ext_fun_expand_cost
        ext_fun_expand_precompute
        ext_fun_expand_dyn
        model_external_shared_lib_dir
        model_external_shared_lib_name

        with_solution_sens_wrt_params
        with_value_sens_wrt_params
        generate_hess
        sens_forw_p
    end

    methods
        function obj = AcadosCodeGenOptions()
            obj.cython_include_dirs = []; % just for python compatibility

            acados_folder = getenv('ACADOS_INSTALL_DIR');
            obj.acados_include_path = [acados_folder, '/include'];
            obj.acados_link_libs = struct();

            if ismac
                obj.shared_lib_ext = '.dylib';
            else
                obj.shared_lib_ext = '.so';
            end

            if ismac
                obj.os = 'mac';
            elseif isunix
                obj.os = 'unix';
            else
                obj.os = 'pc';
            end

            % public
            obj.acados_lib_path = [acados_folder, '/lib'];
            obj.json_file = '';
            obj.code_export_directory = '';
            obj.acados_version = '';
            obj.casadi_code_gen_options = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double', 'force_canonical', false);

            % check whether flags are provided by environment variable
            env_var = getenv("ACADOS_EXT_FUN_COMPILE_FLAGS");
            if isempty(env_var)
                obj.ext_fun_compile_flags = '-O2';
            else
                obj.ext_fun_compile_flags = env_var;
            end
            obj.ext_fun_expand_constr = false;
            obj.ext_fun_expand_cost = false;
            obj.ext_fun_expand_precompute = false;
            obj.ext_fun_expand_dyn = false;

            obj.model_external_shared_lib_dir = [];
            obj.model_external_shared_lib_name = [];

            obj.with_solution_sens_wrt_params = false;
            obj.with_value_sens_wrt_params = false;
            obj.sens_forw_p = false;
            obj.generate_hess = false;
        end

        function make_consistent(obj)

            if ~islogical(obj.sens_forw_p)
                error('sens_forw_p should be a boolean.');
            end

            acados_folder = getenv('ACADOS_INSTALL_DIR');
            addpath(fullfile(acados_folder, 'external', 'jsonlab'));
            libs = loadjson(fileread(fullfile(obj.acados_lib_path, 'link_libs.json')));
            obj.acados_link_libs = orderfields(libs);
            obj.json_file = absolute_path(obj.json_file);

            try % read acados version from git_commit_hash file
                obj.acados_version = fileread(fullfile(obj.acados_lib_path, 'git_commit_hash'));
            catch
                warning('Could not read acados version from git_commit_hash file.');
                obj.acados_version = '';
            end
            if isempty(obj.code_export_directory)
                obj.code_export_directory = 'c_generated_code';
            end
            obj.code_export_directory = absolute_path(obj.code_export_directory);

            if isempty(obj.casadi_code_gen_options)
                obj.casadi_code_gen_options = struct();
            end


            if isfield(obj.casadi_code_gen_options, 'mex') && obj.casadi_code_gen_options.mex
                warning('casadi_code_gen_options.mex is set to true, this is not supported by acados. Setting it to false.');
            end
            if isfield(obj.casadi_code_gen_options, 'casadi_int') && ~strcmp(obj.casadi_code_gen_options.casadi_int, 'int')
                warning('casadi_code_gen_options.casadi_int is set to a value other than "int", this is not supported by acados. Setting it to "int".');
            end
            if isfield(obj.casadi_code_gen_options, 'casadi_real') && ~strcmp(obj.casadi_code_gen_options.casadi_real, 'double')
                warning('casadi_code_gen_options.casadi_real is set to a value other than "double", this is not supported by acados. Setting it to "double".');
            end

            obj.casadi_code_gen_options.mex = false;
            obj.casadi_code_gen_options.casadi_int = 'int';
            obj.casadi_code_gen_options.casadi_real = 'double';
            try
                CodeGenerator('foo', struct('force_canonical', true));
            catch
                if isfield(obj.casadi_code_gen_options, 'force_canonical')
                    warning("CasADi version does not support 'force_canonical' option. Removing it from casadi_code_gen_options.");
                    obj.casadi_code_gen_options = rmfield(obj.casadi_code_gen_options, 'force_canonical');
                end
            end
        end

        function s = to_struct(self)
            if exist('properties')
                publicProperties = eval('properties(self)');
            else
                publicProperties = fieldnames(self);
            end
            s = struct();
            for fi = 1:numel(publicProperties)
                s.(publicProperties{fi}) = self.(publicProperties{fi});
            end

            orderfields(s);
        end
    end
    methods (Static)
        function obj = from_struct(s)
            % Create AcadosCodeGenOptions from a struct (e.g. decoded from JSON).
            obj = AcadosCodeGenOptions();
            fields = fieldnames(s);
            for i = 1:length(fields)
                f = fields{i};
                % direct assignment for simple fields
                try
                    obj.(f) = s.(f);
                catch
                    % ignore unknown fields
                    warning(['Could not assign field ' f ' in AcadosCodeGenOptions.from_struct']);
                end
            end
        end
    end
end
