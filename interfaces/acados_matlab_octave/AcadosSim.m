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

classdef AcadosSim < handle

    properties
        % file structure
        acados_include_path
        acados_lib_path
        shared_lib_ext
        json_file
        cython_include_dirs
        code_export_directory
        % struct / object
        dims
        model
        solver_options
        % plain data
        parameter_values
        problem_class
        external_function_files_model
    end

    methods
        function obj = AcadosSim()
            % most fields are initialized as a placeholder
            obj.acados_include_path = [];
            obj.acados_lib_path = [];
            obj.shared_lib_ext = [];
            obj.json_file = 'acados_sim.json';
            obj.cython_include_dirs = [];
            obj.code_export_directory = 'c_generated_code';

            obj.dims = AcadosSimDims();
            obj.model = AcadosModel();

            obj.solver_options = AcadosSimOptions();

            obj.parameter_values = [];
            obj.problem_class = 'SIM';
        end

        function make_consistent(self)

            % file structures
            acados_folder = getenv('ACADOS_INSTALL_DIR');
            self.acados_include_path = [acados_folder, '/include'];
            self.acados_lib_path = [acados_folder, '/lib'];
            self.shared_lib_ext = '.so';
            if ismac()
                self.shared_lib_ext = '.dylib';
            end
            if isempty(self.json_file)
                self.json_file = 'acados_sim.json';
            end

            % model
            self.model.make_consistent(self.dims);

            if self.dims.np_global > 0
                error('p_global is not supported for AcadosSim.')
            end

            % detect GNSF structure
            if strcmp(self.solver_options.integrator_type, 'GNSF')
                % TODO: interface these options
                gnsf_transcription_opts = struct();
                if self.dims.gnsf_nx1 + self.dims.gnsf_nx2 ~= self.dims.nx
                    detect_gnsf_structure(self.model, self.dims, gnsf_transcription_opts);
                else
                    warning('No GNSF model detected, assuming all required fields are set.')
                end
            end

            % parameters
            if self.dims.np > 0
                if isempty(self.parameter_values)
                    warning(['opts_struct.parameter_values are not set.', ...
                                10 'Using zeros(np,1) by default.' 10 'You can update them later using set().']);
                    self.parameter_values = zeros(self.dims.np,1);
                end
            end

            % solver options checks
            opts = self.solver_options;
            if ~strcmp(opts.integrator_type, 'ERK') && ~strcmp(opts.integrator_type, 'IRK') && ...
               ~strcmp(opts.integrator_type, 'GNSF') && ~strcmp(opts.integrator_type, 'DISCRETE')
                error(['integrator_type = ', opts.integrator_type, ' not available. Choose ERK, IRK, GNSF.']);
            end

            if length(opts.num_stages) ~= 1
                error('num_stages should be a scalar.');
            end
            if length(opts.num_steps) ~= 1
                error('num_steps should be a scalar.');
            end
            if length(opts.newton_iter) ~= 1
                error('newton_iter should be a scalar.');
            end
            % check bool options
            if ~islogical(opts.sens_forw)
                error('sens_forw should be a boolean.');
            end
            if ~islogical(opts.sens_adj)
                error('sens_adj should be a boolean.');
            end
            if ~islogical(opts.sens_algebraic)
                error('sens_algebraic should be a boolean.');
            end
            if ~islogical(opts.sens_hess)
                error('sens_hess should be a boolean.');
            end
            if ~islogical(opts.output_z)
                error('output_z should be a boolean.');
            end
            if ~strcmp(opts.collocation_type, "GAUSS_LEGENDRE") && ~strcmp(opts.collocation_type, "GAUSS_RADAU_IIA")
                error(['collocation_type = ', opts.collocation_type, ' not available. Choose GAUSS_LEGENDRE, GAUSS_RADAU_IIA.']);
            end

            if(strcmp(self.solver_options.integrator_type, "ERK"))
                if(self.solver_options.num_stages == 1 || self.solver_options.num_stages == 2 || ...
                    self.solver_options.num_stages == 3 || self.solver_options.num_stages == 4)
                else
                    error(['ERK: num_stages = ', num2str(self.solver_options.num_stages) ' not available. Only number of stages = {1,2,3,4} implemented!']);
                end
            end
        end

        function generate_external_functions(self)
            if nargin < 2
                % options for code generation
                code_gen_opts = struct();
                code_gen_opts.generate_hess = self.solver_options.sens_hess;
                code_gen_opts.code_export_directory = self.code_export_directory;
                code_gen_opts.ext_fun_expand_dyn = self.solver_options.ext_fun_expand_dyn;
                code_gen_opts.ext_fun_expand_cost = false;
                code_gen_opts.ext_fun_expand_constr = false;
                code_gen_opts.ext_fun_expand_precompute = false;

                context = GenerateContext(self.model.p_global, self.model.name, code_gen_opts);
            else
                code_gen_opts = context.code_gen_opts;
            end

            model_dir = fullfile(pwd, code_gen_opts.code_export_directory, [self.model.name '_model']);
            check_dir_and_create(model_dir);

            if strcmp(self.model.dyn_ext_fun_type, 'generic')
                copyfile(fullfile(pwd, self.model.dyn_generic_source), model_dir);
                context.add_external_function_file(ocp.model.dyn_generic_source, model_dir);

            elseif strcmp(self.model.dyn_ext_fun_type, 'casadi')
                import casadi.*
                check_casadi_version();
                switch self.solver_options.integrator_type
                    case 'ERK'
                        generate_c_code_explicit_ode(context, self.model, model_dir);
                    case 'IRK'
                        generate_c_code_implicit_ode(context, self.model, model_dir);
                    case 'GNSF'
                        generate_c_code_gnsf(context, self.model, model_dir);
                    case 'DISCRETE'
                        error('Discrete dynamics not supported in AcadosSim yet.')
                        % generate_c_code_discrete_dynamics(context, self.model, model_dir);
                    otherwise
                        error('Unknown integrator type.')
                end
            else
                error('Unknown dyn_ext_fun_type.')
            end

            context.finalize();
            self.external_function_files_model = context.get_external_function_file_list(false);
        end

        function dump_to_json(self, json_file)
            if nargin < 2
                json_file = self.json_file;
            end

            %% remove CasADi objects from model
            model.name = self.model.name;
            model.dyn_ext_fun_type = self.model.dyn_ext_fun_type;
            model.dyn_generic_source = self.model.dyn_generic_source;
            model.dyn_disc_fun_jac_hess = self.model.dyn_disc_fun_jac_hess;
            model.dyn_disc_fun_jac = self.model.dyn_disc_fun_jac;
            model.dyn_disc_fun = self.model.dyn_disc_fun;
            model.gnsf_nontrivial_f_LO = self.model.gnsf_nontrivial_f_LO;
            model.gnsf_purely_linear = self.model.gnsf_purely_linear;
            self.model = model;
            % jsonlab
            acados_folder = getenv('ACADOS_INSTALL_DIR');
            addpath(fullfile(acados_folder, 'external', 'jsonlab'))

            %% post process numerical data (mostly cast scalars to 1-dimensional cells)
            % parameter values
            self.parameter_values = reshape(num2cell(self.parameter_values), [1, self.dims.np]);

            %% dump JSON file
            sim_json_struct = self.struct();
            sim_json_struct.dims = self.dims.struct();
            sim_json_struct.solver_options = self.solver_options.struct();

            % add compilation information to json
            libs = loadjson(fileread(fullfile(acados_folder, 'lib', 'link_libs.json')));
            sim_json_struct.acados_link_libs = libs;
            if ismac
                sim_json_struct.os = 'mac';
            elseif isunix
                sim_json_struct.os = 'unix';
            else
                sim_json_struct.os = 'pc';
            end

            json_string = savejson('', sim_json_struct, 'ForceRootName', 0);

            fid = fopen(self.json_file, 'w');
            if fid == -1, error('Cannot create JSON file'); end
            fwrite(fid, json_string, 'char');
            fclose(fid);
        end

        function render_templates(self)

            json_fullfile = fullfile(pwd, self.json_file);

            acados_root_dir = getenv('ACADOS_INSTALL_DIR');
            acados_template_folder = fullfile(acados_root_dir,...
                                'interfaces', 'acados_template', 'acados_template');

            t_renderer_location = get_tera();

            %% load json data
            acados_sim = loadjson(fileread(json_fullfile));
            model_name = acados_sim.model.name;

            %% render templates
            matlab_template_path = 'matlab_templates';
            main_dir = pwd;
            chdir(self.code_export_directory);

            % cell array with entries (template_file, output file)
            template_list = { ...
                {'main_sim.in.c', ['main_sim_', model_name, '.c']}, ...
                {fullfile(matlab_template_path, 'mex_sim_solver.in.m'), [model_name, '_mex_sim_solver.m']}, ...
                {fullfile(matlab_template_path, 'make_mex_sim.in.m'), ['make_mex_sim_', model_name, '.m']}, ...
                {fullfile(matlab_template_path, 'acados_sim_create.in.c'), ['acados_sim_create_', model_name, '.c']}, ...
                {fullfile(matlab_template_path, 'acados_sim_free.in.c'), ['acados_sim_free_', model_name, '.c']}, ...
                {fullfile(matlab_template_path, 'acados_sim_set.in.c'), ['acados_sim_set_', model_name, '.c']}, ...
                {'acados_sim_solver.in.c', ['acados_sim_solver_', model_name, '.c']}, ...
                {'acados_sim_solver.in.h', ['acados_sim_solver_', model_name, '.h']}, ...
                {fullfile(matlab_template_path, 'acados_sim_solver_sfun.in.c'), ['acados_sim_solver_sfunction_', model_name, '.c']}, ...
                {fullfile(matlab_template_path, 'make_sfun_sim.in.m'), ['make_sfun_sim_', model_name, '.m']}, ...
                {'Makefile.in', 'Makefile'}, ...
                {'CMakeLists.in.txt', 'CMakeLists.txt'}};

            num_entries = length(template_list);
            for n=1:num_entries
                entry = template_list{n};
                render_file( entry{1}, entry{2}, json_fullfile);
            end

            c_dir = pwd;
            chdir([model_name, '_model']);
            render_file( 'model.in.h', [model_name, '_model.h'], json_fullfile);
            cd(c_dir);

            fprintf('Successfully rendered acados templates!\n');
            cd(main_dir)
        end

        function s = struct(self)
            if exist('properties')
                publicProperties = eval('properties(self)');
            else
                publicProperties = fieldnames(self);
            end
            s = struct();
            for fi = 1:numel(publicProperties)
                s.(publicProperties{fi}) = self.(publicProperties{fi});
            end
        end
    end

end % class
