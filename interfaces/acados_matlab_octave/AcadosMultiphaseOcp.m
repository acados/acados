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

classdef AcadosMultiphaseOcp < handle
    properties
        N_list
        n_phases
        N_horizon

        phases_dims
        cost
        constraints
        solver_options
        mocp_opts
        dummy_ocp_list

        model
        parameter_values % initial value of the parameter
        p_global_values % initial value of global parameters
        problem_class
        simulink_opts
        name

        cython_include_dirs
        code_export_directory
        json_file
        shared_lib_ext
        acados_include_path
        acados_lib_path

        % detected fields
        start_idx
        end_idx
        cost_start_idx

        external_function_files_ocp
        external_function_files_model
    end
    methods
        function obj = AcadosMultiphaseOcp(N_list)
            if any(N_list < 1)
                error('N_list must be a list of positive integers');
            end
            n_phases = length(N_list);
            obj.n_phases = n_phases;
            obj.N_list = N_list;
            obj.N_horizon = sum(N_list);

            obj.phases_dims = cell(n_phases, 1);
            obj.cost = cell(n_phases, 1);
            obj.constraints = cell(n_phases, 1);
            obj.model = cell(n_phases, 1);

            for i=1:n_phases
                obj.phases_dims{i} = AcadosOcpDims();
                obj.cost{i} = AcadosOcpCost();
                obj.constraints{i} = AcadosOcpConstraints();
                obj.model{i} = AcadosModel();
            end

            obj.dummy_ocp_list = cell(n_phases, 1);

            obj.solver_options = AcadosOcpOptions();
            obj.solver_options.N_horizon = obj.N_horizon; % NOTE: to not change options when making ocp consistent

            obj.mocp_opts = AcadosMultiphaseOptions();

            obj.parameter_values = cell(n_phases, 1);
            obj.p_global_values = [];
            obj.problem_class = 'MOCP';
            obj.simulink_opts = [];
            obj.cython_include_dirs = [];
            obj.json_file = 'mocp.json';
            obj.shared_lib_ext = '.so';
            obj.name = 'ocp';
            if ismac()
                obj.shared_lib_ext = '.dylib';
            end
            obj.code_export_directory = 'c_generated_code';

            % set include and lib path
            acados_folder = getenv('ACADOS_INSTALL_DIR');
            obj.acados_include_path = [acados_folder, '/include'];
            obj.acados_lib_path = [acados_folder, '/lib'];

        end


        function set_phase(self, ocp, phase_idx)
            % Note: phase_idx is 1-indexed in contrast to Python!
            if phase_idx > self.n_phases
                error('phase_idx must be less than or equal to the number of phases');
            end

            % Check solver options
            non_default_opts = find_non_default_fields_of_obj(ocp.solver_options);

            if ~isempty(non_default_opts)
                fprintf('WARNING: set_phase: Phase %d contains non-default solver options: %s, which will be ignored.\n', ...
                        phase_idx, strjoin(non_default_opts, ', '));
                fprintf('Solver options need to be set via AcadosMultiphaseOcp.solver_options and via AcadosMultiphaseOcp.mocp_opts for options that can only vary in MOCP.\n');
            end

            % set phase
            self.model{phase_idx} = ocp.model;
            self.cost{phase_idx} = ocp.cost;
            self.constraints{phase_idx} = ocp.constraints;
            self.parameter_values{phase_idx} = ocp.parameter_values;

            if ~isempty(self.p_global_values)
                fprintf('WARNING: set_phase: Phase %d contains p_global values which will be ignored.\n', phase_idx);
            end
        end

        function make_consistent(self)
            % check options
            self.mocp_opts.make_consistent(self.solver_options, self.n_phases);

            % check phases formulation objects are distinct
            if ~is_octave() % octave does not support object comparison
                for i=1:self.n_phases
                    for j=i+1:self.n_phases
                        if self.model{i} == self.model{j}
                            error('model objects must be distinct for each phase');
                        end
                        if self.cost{i} == self.cost{j}
                            error('cost objects must be distinct for each phase');
                        end
                        if self.constraints{i} == self.constraints{j}
                            error('constraints objects must be distinct for each phase');
                        end
                    end
                end
            end

            % check N_horizon
            if self.N_horizon ~= sum(self.N_list)
                error('N_horizon must be equal to the sum of N_list, N_horizon is detected automatically for AcadosMultiphaseOcp and should not be set manually.');
            end

            % compute phase indices
            phase_idx = cumsum([0, self.N_list]);
            self.start_idx = phase_idx(1:end-1);
            self.end_idx = phase_idx(2:end);

            self.cost_start_idx = phase_idx;
            self.cost_start_idx(1) = self.cost_start_idx(1) + 1;

            % make model names unique if necessary
            model_name_list = cell(self.n_phases, 1);
            for i=1:self.n_phases
                model_name_list{i} = self.model{i}.name;
            end
            n_names = length(unique(model_name_list));
            if n_names ~= self.n_phases
                disp('model names are not unique: got');
                disp(model_name_list);
                disp('adding _i to model names');
                for i=1:self.n_phases
                    self.model{i}.name = [self.model{i}.name, '_', num2str(i)];
                end
                model_name_list = cell(self.n_phases, 1);
                for i=1:self.n_phases
                    model_name_list{i} = self.model{i}.name;
                end
                disp('new model names are');
                disp(model_name_list);
            end

            % make phase OCPs consistent, warn about unused fields
            for i=1:self.n_phases
                ocp = AcadosOcp();
                ocp.dims = self.phases_dims{i};
                ocp.model = self.model{i};
                ocp.constraints = self.constraints{i};
                ocp.cost = self.cost{i};
                ocp.parameter_values = self.parameter_values{i};
                ocp.p_global_values = self.p_global_values;
                ocp.solver_options = self.solver_options;

                % set phase dependent options
                ocp.solver_options.integrator_type = self.mocp_opts.integrator_type{i};
                ocp.solver_options.collocation_type = self.mocp_opts.collocation_type{i};
                ocp.solver_options.cost_discretization = self.mocp_opts.cost_discretization{i};

                % check for non-default fields in terminal/initial phase that are not used
                if i ~= self.n_phases % not terminal phase
                    nondefault_fields = {};

                    nondefault_fields = [nondefault_fields, find_non_default_fields_of_obj(ocp.cost, 'terminal')];
                    nondefault_fields = [nondefault_fields, find_non_default_fields_of_obj(ocp.constraints, 'terminal')];
                    nondefault_fields = [nondefault_fields, find_non_default_fields_of_obj(ocp.model, 'terminal')];

                    if ~isempty(nondefault_fields)
                        disp(['Phase ', num2str(i), ' contains non-default terminal fields: ', strjoin(nondefault_fields, ', '), ', which will be ignored.']);
                    end
                elseif i ~= 1 % not initial phase
                    nondefault_fields = {};

                    nondefault_fields = [nondefault_fields, find_non_default_fields_of_obj(ocp.cost, 'initial')];
                    nondefault_fields = [nondefault_fields, find_non_default_fields_of_obj(ocp.constraints, 'initial')];
                    nondefault_fields = [nondefault_fields, find_non_default_fields_of_obj(ocp.model, 'initial')];

                    if ~isempty(nondefault_fields)
                        disp(['Phase ', num2str(i), ' contains non-default initial fields: ', strjoin(nondefault_fields, ', '), ', which will be ignored.']);
                    end
                end

                disp(['Calling make_consistent for phase ', num2str(i), '.']);
                ocp.make_consistent(true);
                % use the updated objects that are not handles
                self.parameter_values{i} = ocp.parameter_values;

                self.dummy_ocp_list{i} = ocp;
            end


            % check for transition consistency
            nx_list = zeros(self.n_phases, 1);
            for i=1:self.n_phases
                nx_list(i) = self.phases_dims{i}.nx;
            end
            for i=2:self.n_phases
                if nx_list(i) ~= nx_list(i-1)
                    if self.phases_dims{i-1}.nx_next ~= self.phases_dims{i}.nx
                        error(['detected stage transition with different nx from phase ', num2str(i-1), ' to ', num2str(i), ', which is only supported for nx_next = nx, got nx_next = ', num2str(self.phases_dims{i-1}.nx_next), ' and nx = ', num2str(self.phases_dims{i}.nx), '.']);
                    end
                    if self.N_list(i-1) ~= 1 || ~strcmp(self.mocp_opts.integrator_type{i-1}, 'DISCRETE')
                        error(['detected stage transition with different nx from phase ', num2str(i-1), ' to ', num2, ', which is only supported for integrator_type=''DISCRETE'' and N_list[i] == 1.']);
                    end
                end
            end
        end

        function template_list = get_template_list(self)
            % returns a cell of cells in the form:
            % (input_filename, output_filname)
            % or
            % (input_filename, output_filname, output_directory)

            template_list = {};
            template_list{end+1} = {'main_multi.in.c', ['main_', self.name, '.c']};
            template_list{end+1} = {'acados_multi_solver.in.h', ['acados_solver_', self.name, '.h']};
            template_list{end+1} = {'acados_multi_solver.in.c', ['acados_solver_', self.name, '.c']};
            template_list{end+1} = {'multi_CMakeLists.in.txt', 'CMakeLists.txt'};
            template_list{end+1} = {'multi_Makefile.in', 'Makefile'};

            % MEX files
            matlab_template_path = 'matlab_templates';
            template_list{end+1} = {fullfile(matlab_template_path, 'mex_solver.in.m'), [self.name, '_mex_solver.m']};
            template_list{end+1} = {fullfile(matlab_template_path, 'make_mex.in.m'), ['make_mex_', self.name, '.m']};
            template_list{end+1} = {fullfile(matlab_template_path, 'acados_mex_create.in.c'), ['acados_mex_create_', self.name, '.c']};
            template_list{end+1} = {fullfile(matlab_template_path, 'acados_mex_free.in.c'), ['acados_mex_free_', self.name, '.c']};
            template_list{end+1} = {fullfile(matlab_template_path, 'acados_mex_solve.in.c'), ['acados_mex_solve_', self.name, '.c']};
            template_list{end+1} = {fullfile(matlab_template_path, 'acados_mex_set.in.c'), ['acados_mex_set_', self.name, '.c']};
            if self.phases_dims{1}.n_global_data > 0
                template_list{end+1} = {'p_global_precompute_fun.in.h',  [self.name, '_p_global_precompute_fun.h']};
            end
            % Simulink
            if ~isempty(self.simulink_opts)
                template_list{end+1} = {fullfile(matlab_template_path, 'acados_solver_sfun.in.c'), ['acados_solver_sfunction_', self.name, '.c']};
                template_list{end+1} = {fullfile(matlab_template_path, 'make_sfun.in.m'), ['make_sfun.m']};
                % TODO: do we want to generate simulink sfun for sim solver?
                % if ~strcmp(self.solver_options.integrator_type, 'DISCRETE')
                %     template_list{end+1} = {fullfile(matlab_template_path, 'acados_sim_solver_sfun.in.c'), ['acados_sim_solver_sfunction_', self.name, '.c']};
                %     template_list{end+1} = {fullfile(matlab_template_path, 'make_sfun_sim.in.m'), ['make_sfun_sim.m']};
                % end
                if self.simulink_opts.inputs.rti_phase && self.solver_options.nlp_solver_type ~= 'SQP_RTI'
                    error('rti_phase is only supported for SQP_RTI');
                end
                inputs = self.simulink_opts.inputs;
                nonsupported_mocp_inputs = {'y_ref', 'lg', 'ug', 'cost_W_0', 'cost_W', 'cost_W_e'};
                for i=1:length(nonsupported_mocp_inputs)
                    if inputs.(nonsupported_mocp_inputs{i})
                        error(['Simulink inputs ', nonsupported_mocp_inputs{i}, ' are not supported for MOCP.']);
                    end
                end
            else
                disp("not rendering Simulink related templates, as simulink_opts are not specified.")
            end
        end

        function context = generate_external_functions(self)
            % generate external functions
            code_gen_opts = struct();
            code_gen_opts.generate_hess = strcmp(self.solver_options.hessian_approx, 'EXACT');
            code_gen_opts.with_solution_sens_wrt_params = self.solver_options.with_solution_sens_wrt_params;
            code_gen_opts.with_value_sens_wrt_params = self.solver_options.with_value_sens_wrt_params;
            code_gen_opts.code_export_directory = self.code_export_directory;

            code_gen_opts.ext_fun_expand_dyn = self.solver_options.ext_fun_expand_dyn;
            code_gen_opts.ext_fun_expand_cost = self.solver_options.ext_fun_expand_cost;
            code_gen_opts.ext_fun_expand_constr = self.solver_options.ext_fun_expand_constr;
            code_gen_opts.ext_fun_expand_precompute = self.solver_options.ext_fun_expand_precompute;
            context = GenerateContext(self.model{1}.p_global, self.name, code_gen_opts);

            for i=1:self.n_phases
                disp(['generating external functions for phase ', num2str(i)]);
                if i ~= self.n_phases
                    ignore_terminal = true;
                else
                    ignore_terminal = false;
                end

                if i ~= 1
                    ignore_initial = true;
                else
                    ignore_initial = false;
                end

                % this is the only option that can vary and influence external functions to be generated
                self.dummy_ocp_list{i}.solver_options.integrator_type = self.mocp_opts.integrator_type{i};
                self.dummy_ocp_list{i}.code_export_directory = self.code_export_directory;
                context = self.dummy_ocp_list{i}.setup_code_generation_context(context, ignore_initial, ignore_terminal);
            end

            context.finalize();
            self.external_function_files_model = context.get_external_function_file_list(false);
            self.external_function_files_ocp = context.get_external_function_file_list(true);

            for i=1:self.n_phases
                self.phases_dims{i}.n_global_data = context.get_n_global_data();
            end
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
            % delete keys that should not be used
            s = rmfield(s, 'dummy_ocp_list');
            s.solver_options = self.solver_options.struct();
            s.solver_options = rmfield(s.solver_options, 'integrator_type');
            s.solver_options = rmfield(s.solver_options, 'collocation_type');
            s.solver_options = rmfield(s.solver_options, 'cost_discretization');
        end

        function dump_to_json(self)
            out_struct = orderfields(self.struct());

            % add compilation information to json
            acados_folder = getenv('ACADOS_INSTALL_DIR');
            libs = loadjson(fileread(fullfile(acados_folder, 'lib', 'link_libs.json')));
            out_struct.acados_link_libs = orderfields(libs);
            if ismac
                out_struct.os = 'mac';
            elseif isunix
                out_struct.os = 'unix';
            else
                out_struct.os = 'pc';
            end

            % prepare struct for json dump
            out_struct.p_global_values = reshape(num2cell(self.p_global_values), [1, self.phases_dims{1}.np_global]);
            for i=1:self.n_phases
                out_struct.parameter_values{i} = reshape(num2cell(self.parameter_values{i}), [1, self.phases_dims{i}.np]);
                out_struct.model{i} = orderfields(self.model{i}.convert_to_struct_for_json_dump());
                out_struct.phases_dims{i} = orderfields(self.phases_dims{i}.struct());
                out_struct.cost{i} = orderfields(self.cost{i}.convert_to_struct_for_json_dump());
                out_struct.constraints{i} = orderfields(self.constraints{i}.convert_to_struct_for_json_dump());
            end
            out_struct.solver_options = orderfields(self.solver_options.convert_to_struct_for_json_dump(self.N_horizon));
            out_struct.mocp_opts = orderfields(self.mocp_opts.struct());

            vector_fields = {'model', 'phases_dims', 'cost', 'constraints', 'parameter_values', 'p_global_values'};
            out_struct = prepare_struct_for_json_dump(out_struct, vector_fields, {});

            % add full path to json file
            self.json_file = fullfile(pwd, self.json_file);
            % actual json dump
            json_string = savejson('', out_struct, 'ForceRootName', 0);
            fid = fopen(self.json_file, 'w');
            if fid == -1, error('Cannot create JSON file'); end
            fwrite(fid, json_string, 'char');
            fclose(fid);
        end

        function render_templates(self)

            main_dir = pwd;
            chdir(self.code_export_directory);

            % model templates
            for i=1:self.n_phases
                % this is the only option that can vary and influence external functions to be generated
                self.dummy_ocp_list{i}.solver_options.integrator_type = self.mocp_opts.integrator_type{i};

                template_list = self.dummy_ocp_list{i}.get_external_function_header_templates();
                % dump dummy_ocp
                tmp_json_file = 'tmp_ocp.json';
                self.dummy_ocp_list{i}.dump_to_json(tmp_json_file);
                tmp_json_path = fullfile(pwd, tmp_json_file);

                for j = 1:length(template_list)
                    in_file = template_list{j}{1};
                    out_file = template_list{j}{2};
                    if length(template_list{j}) == 3
                        out_dir = template_list{j}{3};
                        if ~(exist(out_dir, 'dir'))
                            mkdir(out_dir);
                        end
                        out_file = fullfile(out_dir, out_file);
                    end
                    render_file( in_file, out_file, tmp_json_path );
                end
            end
            disp('rendered model templates successfully');

            % check json file
            if ~(exist(self.json_file, 'file'))
                error(['Path "', self.json_file, '" not found!']);
            end

            % solver templates
            template_list = self.get_template_list();

            % Render templates
            for i = 1:length(template_list)
                in_file = template_list{i}{1};
                out_file = template_list{i}{2};
                if length(template_list{i}) == 3
                    out_dir = template_list{i}{3};
                    if ~(exist(out_dir, 'dir'))
                        mkdir(out_dir);
                    end
                    out_file = fullfile(out_dir, out_file);
                end
                render_file( in_file, out_file, self.json_file );
            end

            disp('rendered solver templates successfully!');
            cd(main_dir);
        end
    end % methods
end

