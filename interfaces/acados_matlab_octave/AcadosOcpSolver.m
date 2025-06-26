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

classdef AcadosOcpSolver < handle

    properties (Access = public)
        ocp % MATLAB class AcadosOcp describing the OCP formulation
    end % properties

    properties (Access = private)
        fields = {'x', 'u', 'z', 'sl', 'su', 'lam', 'pi'};
        qp_gettable_fields = {'qp_Q', 'qp_R', 'qp_S', 'qp_q', 'qp_r', 'qp_A', 'qp_B', 'qp_b', 'qp_C', 'qp_D', 'qp_lg', 'qp_ug', 'qp_lbx', 'qp_ubx', 'qp_lbu', 'qp_ubu', 'qp_zl', 'qp_zu', 'qp_Zl', 'qp_Zu'}
        t_ocp % templated solver

        % required info loaded from json
        N_horizon
        solver_options
        problem_class
        name
        has_x0
        nsbu_0
        nbxe_0
    end
    methods

        function obj = AcadosOcpSolver(ocp, varargin)
            %% optional arguments:
            % varargin{1}: solver_creation_opts: this is a struct in which some of the fields can be defined to overwrite the default values.
            % The fields are:
            % - json_file: path to the json file containing the ocp description
            % - build: boolean, if true, the problem specific shared library is compiled
            % - generate: boolean, if true, the C code is generated
            % - compile_mex_wrapper: boolean, if true, the mex wrapper is compiled
            % - compile_interface: can be [], true or false. If [], the interface is compiled if it does not exist.
            % - output_dir: path to the directory where the MEX interface is compiled
            obj.ocp = ocp;

            % optional arguments
            % solver creation options
            default_solver_creation_opts = struct('json_file', '', ...
                    'build', true, ...
                    'generate', true, ...
                    'compile_mex_wrapper', true, ...
                    'compile_interface', [], ...
                    'output_dir', fullfile(pwd, 'build'));
            if length(varargin) > 0
                solver_creation_opts = varargin{1};
                % set non-specified opts to default
                fields = fieldnames(default_solver_creation_opts);
                for i = 1:length(fields)
                    if ~isfield(solver_creation_opts, fields{i})
                        solver_creation_opts.(fields{i}) = default_solver_creation_opts.(fields{i});
                    end
                end
            else
                solver_creation_opts = default_solver_creation_opts;
            end

            if isempty(ocp) && isempty(solver_creation_opts.json_file)
                error('AcadosOcpSolver: provide either an OCP object or a json file');
            end

            if isempty(ocp)
                json_file = solver_creation_opts.json_file;
            else
                % formulation provided
                if ~isempty(solver_creation_opts.json_file)
                    ocp.json_file = solver_creation_opts.json_file;
                end
                json_file = ocp.json_file;
                if ~isempty(ocp.solver_options.compile_interface) && ~isempty(solver_creation_opts.compile_interface)
                    error('AcadosOcpSolver: provide either compile_interface in OCP object or solver_creation_opts');
                end
                if ~isempty(ocp.solver_options.compile_interface)
                    solver_creation_opts.compile_interface = ocp.solver_options.compile_interface;
                end
                % make consistent
                ocp.make_consistent();
            end

            %% compile mex interface if needed
            obj.compile_mex_interface_if_needed(solver_creation_opts);

            %% generate
            if solver_creation_opts.generate
                obj.generate();
            end

            %% load json, store options in object
            acados_ocp_struct = loadjson(fileread(json_file), 'SimplifyCell', 0);
            obj.problem_class = acados_ocp_struct.problem_class;
            obj.solver_options = acados_ocp_struct.solver_options;
            obj.N_horizon = acados_ocp_struct.solver_options.N_horizon;
            obj.name = acados_ocp_struct.name;

            if strcmp(obj.problem_class, "OCP")
                obj.has_x0 = acados_ocp_struct.constraints.has_x0;
                obj.nsbu_0 = acados_ocp_struct.dims.nsbu;
                obj.nbxe_0 = acados_ocp_struct.dims.nbxe_0;
            elseif strcmp(obj.problem_class, "MOCP")
                obj.has_x0 = acados_ocp_struct.constraints{1}.has_x0;
                obj.nsbu_0 = acados_ocp_struct.phases_dims{1}.nsbu;
                obj.nbxe_0 = acados_ocp_struct.phases_dims{1}.nbxe_0;
            end
            code_export_directory = acados_ocp_struct.code_export_directory;

            %% compile problem specific shared library
            if solver_creation_opts.build
                obj.compile_ocp_shared_lib(code_export_directory);
            end

            %% create solver
            return_dir = pwd();
            cd(code_export_directory)

            mex_solver_name = sprintf('%s_mex_solver', obj.name);
            mex_solver = str2func(mex_solver_name);
            obj.t_ocp = mex_solver(solver_creation_opts);
            addpath(pwd());

            cd(return_dir);
        end

        function solve(obj)
            obj.t_ocp.solve();
        end

        % TODO: remove this! in v.0.5.0!
        function generate_c_code(obj, simulink_opts)
            warning('acados_ocp will be deprecated in the future. Use AcadosOcpSolver instead. For more information on the major acados MATLAB interface overhaul, see https://github.com/acados/acados/releases/tag/v0.4.0');
            if nargin < 2
                warning("Code is generated in the constructor of AcadosOcpSolver.")
            else
                error("If you want to provide simulink options, provide them in AcadosOcp.simulink_opts.")
            end
        end

        function eval_param_sens(obj, field, stage, index)
            obj.t_ocp.eval_param_sens(field, stage, index);
        end

        function value = get_cost(obj)
            value = obj.t_ocp.get_cost();
        end

        function value = evaluate_constraints_and_get_violation(obj)
            % returns the constraint violations for all stages in a cell array.
            % values > 0 indicate a violation of the constraints.
            value = obj.t_ocp.evaluate_constraints_and_get_violation();
        end

        function violation_idx = get_constraint_indices_with_violation(obj, tol)
            % computes the indices of the constraints that are violated with respect to the given tolerance, in the form
            % `[stage_index_0, constraint_index_0,
            % stage_index_1, constraint_index_1,
            % ....]
            % all indices are zero -based.
            ineq_fun = obj.evaluate_constraints_and_get_violation();
            if nargin < 2
                tol = obj.solver_options.nlp_solver_tol_ineq;
            end
            violation_idx = [];
            for i=1:length(ineq_fun)
                idx_i = find(ineq_fun{i} >= tol);
                violation_idx = [violation_idx; [(i-1)*ones(length(idx_i), 1), idx_i-1]];
            end
        end

        function set(obj, field, value, varargin)
            obj.t_ocp.set(field, value, varargin{:});
        end

        function status = custom_update(obj, data)
            status = obj.t_ocp.custom_update(data);
        end

        function value = get_qp_scaling_constraints(obj, stage)
            % returns the qp scaling constraints for the given stage
            value = obj.t_ocp.get('qpscaling_constr', stage);
        end

        function value = get_qp_scaling_objective(obj)
            % returns the qp scaling objective
            value = obj.t_ocp.get('qpscaling_obj');
        end

        function value = get(obj, field, varargin)

            if strcmp('hess_block', field)

                if length(varargin) > 0
                    n = varargin{1};

                    if n < obj.solver_options.N_horizon
                        Q = obj.get('qp_Q', n);
                        R = obj.get('qp_R', n);
                        S = obj.get('qp_S', n);

                        value = [R, S; S', Q];
                    else
                        value = obj.get('qp_Q', n);
                    end

                    return;
                else
                    value = cell(obj.solver_options.N_horizon, 1);
                    for n=0:(obj.solver_options.N_horizon-1)
                        Q = obj.get('qp_Q', n);
                        R = obj.get('qp_R', n);
                        S = obj.get('qp_S', n);

                        value{n+1} = [R, S; S', Q];
                    end
                    value{end+1} = obj.get('qp_Q', obj.solver_options.N_horizon);
                    return;
                end
            elseif strcmp('pc_hess_block', field)

                if ~strncmp(obj.solver_options.qp_solver, 'PARTIAL_CONDENSING', length('PARTIAL_CONDENSING'))
                    error("Getting hessian block of partially condensed QP only works for PARTIAL_CONDENSING QP solvers");
                end
                if length(varargin) > 0
                    n = varargin{1};
                    all_blocks = obj.get('qp_solver_cond_H');
                    value = all_blocks{n+1};
                    return;
                else
                    value = obj.get('qp_solver_cond_H');
                    return;
                end
            else
                value = obj.t_ocp.get(field, varargin{:});
            end

            % make symmetric (only lower triangular stored internally)
            if strcmp('qp_Q', field) || strcmp('qp_R', field)
                if iscell(value)
                    for i=1:length(value)
                        if length(value{i}) > 1
                            value{i} = tril(value{i}) + tril(value{i}, -1)';
                        end
                    end
                else
                    if length(value) > 1
                        value = tril(value) + tril(value, -1)';
                    end
                end
            end
        end

        function [] = store_iterate(obj, varargin)
            %%%  Stores the current iterate of the ocp solver in a json file.
            %%% param1: filename: if not set, use model_name + timestamp + '.json'
            %%% param2: overwrite: if false and filename exists add timestamp to filename

            obj.t_ocp.store_iterate(varargin{:});
        end


        function [] = load_iterate(obj, filename)
            obj.t_ocp.load_iterate(filename);
        end

        function [] = load_iterate_from_obj(obj, iterate)
            %%%  Loads the iterate from an AcadosOcpIterate object.
            %%% param1: iterate: AcadosOcpIterate object containing the iterate to load

            if ~isa(iterate, 'AcadosOcpIterate')
                error('load_iterate_from_obj: iterate needs to be of type AcadosOcpIterate');
            end

            for i = 1:obj.solver_options.N_horizon + 1
                obj.t_ocp.set('x', iterate.x_traj{i, 1}, i-1);
                obj.t_ocp.set('sl', iterate.sl_traj{i, 1}, i-1);
                obj.t_ocp.set('su', iterate.su_traj{i, 1}, i-1);
                obj.t_ocp.set('lam', iterate.lam_traj{i, 1}, i-1);
            end
            for i = 1:obj.solver_options.N_horizon
                obj.t_ocp.set('u', iterate.u_traj{i, 1}, i-1);
                obj.t_ocp.set('pi', iterate.pi_traj{i, 1}, i-1);
                if ~isempty(iterate.z_traj{i, 1})
                    obj.t_ocp.set('z', iterate.z_traj{i, 1}, i-1);
                end
            end
        end
        function iterate = get_iterate(obj, iteration)
            if iteration > obj.get('nlp_iter')
                error("iteration needs to be nonnegative and <= nlp_iter.");
            end

            if ~obj.solver_options.store_iterates
                error("get_iterate: the solver option store_iterates needs to be true in order to get iterates.");
            end

            if strcmp(obj.solver_options.nlp_solver_type, 'SQP_RTI')
                error("get_iterate: SQP_RTI not supported.");
            end

            N_horizon = obj.solver_options.N_horizon;

            x_traj = cell(N_horizon + 1, 1);
            u_traj = cell(N_horizon, 1);
            z_traj = cell(N_horizon, 1);
            sl_traj = cell(N_horizon + 1, 1);
            su_traj = cell(N_horizon + 1, 1);
            pi_traj = cell(N_horizon, 1);
            lam_traj = cell(N_horizon + 1, 1);

            for n=1:N_horizon
                x_traj{n, 1} = obj.t_ocp.get('x', n-1, iteration);
                u_traj{n, 1} = obj.t_ocp.get('u', n-1, iteration);
                z_traj{n, 1} = obj.t_ocp.get('z', n-1, iteration);
                sl_traj{n, 1} = obj.t_ocp.get('sl', n-1, iteration);
                su_traj{n, 1} = obj.t_ocp.get('su', n-1, iteration);
                pi_traj{n, 1} = obj.t_ocp.get('pi', n-1, iteration);
                lam_traj{n, 1} = obj.t_ocp.get('lam', n-1, iteration);
            end

            x_traj{N_horizon+1, 1} = obj.t_ocp.get('x', N_horizon, iteration);
            sl_traj{N_horizon+1, 1} = obj.t_ocp.get('sl', N_horizon, iteration);
            su_traj{N_horizon+1, 1} = obj.t_ocp.get('su', N_horizon, iteration);
            lam_traj{N_horizon+1, 1} = obj.t_ocp.get('lam', N_horizon, iteration);

            iterate = AcadosOcpIterate(x_traj, u_traj, z_traj, ...
                    sl_traj, su_traj, pi_traj, lam_traj);
        end

        function iterates = get_iterates(obj)
            nlp_iter = obj.get('nlp_iter');
            iterates_cell = cell(nlp_iter+1, 1);

            for n=1:(nlp_iter+1)
                iterates_cell{n} = obj.get_iterate(n-1);
            end

            iterates = AcadosOcpIterates(iterates_cell);
        end

        function print(obj, varargin)
            obj.t_ocp.print(varargin{:});
        end


        function reset(obj)
            obj.t_ocp.reset();
        end


        function result = qp_diagnostics(obj, varargin)
            % Compute some diagnostic values for the last QP.
            % result = ocp_solver.qp_diagnostics([partially_condensed_qp=false])

            % returns a struct with the following fields:
            % - min_eigv_stage: array with minimum eigenvalue for each Hessian block.
            % - max_eigv_stage: array with maximum eigenvalue for each Hessian block.
            % - condition_number_stage: array with condition number for each Hessian block.
            % - min_eigv_global: minimum eigenvalue for the full Hessian.
            % - min_abs_eigv_global: minimum absolute eigenvalue for the full Hessian.
            % - max_eigv_global: maximum eigenvalue for the full Hessian.
            % - max_abs_eigv_global: maximum absolute eigenvalue for the full Hessian.
            % - condition_number_global: condition number for the full Hessian.

            if length(varargin) > 0
                partially_condensed_qp = varargin{1};
            else
                partially_condensed_qp = false;
            end

            if partially_condensed_qp
                num_blocks = obj.solver_options.qp_solver_cond_N + 1;
            else
                num_blocks = obj.N_horizon + 1;
            end
            result = struct();
            result.min_eigv_stage = zeros(num_blocks, 1);
            result.max_eigv_stage = zeros(num_blocks, 1);
            result.condition_number_stage = zeros(num_blocks, 1);
            min_abs_val = inf;
            max_abs_val = -inf;
            max_ev = -inf;
            min_ev = inf;

            for n=1:num_blocks
                if partially_condensed_qp
                    hess_block = obj.get('pc_hess_block', n-1);
                else
                    hess_block = obj.get('hess_block', n-1);
                end
                eigvals = eig(hess_block);
                max_ev = max(max_ev, max(eigvals));
                max_abs_val = max(max_abs_val, max(abs(eigvals)));
                min_ev = min(min_ev, min(eigvals));
                min_abs_val = min(min_abs_val, min(abs(eigvals)));
                result.min_eigv_stage(n) = min(eigvals);
                result.max_eigv_stage(n) = max(eigvals);
                result.condition_number_stage(n) = max(eigvals) / min(eigvals);
            end
            result.condition_number_global = max_abs_val / min_abs_val;
            result.max_eigv_global = max_ev;
            result.max_abs_eigv_global = max_abs_val;
            result.min_eigv_global = min_ev;
            result.min_abs_eigv_global = min_abs_val;
        end


        function dump_last_qp_to_json(obj, filename)
            qp_data = struct();

            lN = length(num2str(obj.solver_options.N_horizon+1));
            n_fields = length(obj.qp_gettable_fields);
            for n=1:n_fields

                field = obj.qp_gettable_fields{n};
                for i=0:obj.solver_options.N_horizon-1
                    s_indx = sprintf(strcat('%0', num2str(lN), 'd'), i);
                    key = strcat(field, '_', s_indx);
                    val = obj.get(field, i);
                    qp_data = setfield(qp_data, key, val);
                end
                if strcmp(field, 'qp_Q') || ...
                   strcmp(field, 'qp_q') || ...
                   strcmp(field, 'qp_C') || ...
                   strcmp(field, 'qp_lg') || ...
                   strcmp(field, 'qp_ug') || ...
                   strcmp(field, 'qp_lbx') || ...
                   strcmp(field, 'qp_ubx') || ...
                   strcmp(field, 'qp_zl') || ...
                   strcmp(field, 'qp_zu') || ...
                   strcmp(field, 'qp_Zl') || ...
                   strcmp(field, 'qp_Zu')
                    s_indx = sprintf(strcat('%0', num2str(lN), 'd'), obj.solver_options.N_horizon);
                    key = strcat(field, '_', s_indx);
                    val = obj.get(field, obj.solver_options.N_horizon);
                    qp_data = setfield(qp_data, key, val);
                end
            end

            % save
            json_string = savejson('', qp_data, 'ForceRootName', 0, struct('FloatFormat', '%.5f'));

            fid = fopen(filename, 'w');
            if fid == -1, error('Cannot create json file'); end
            fwrite(fid, json_string, 'char');
            fclose(fid);
        end


        function set_params_sparse(obj, varargin)
            % usage:
            % ocp.set_params_sparse(idx_values, param_values, Optional[stage])
            % updates the parameters with indices idx_values (0 based) at stage with the new values new_p_values.
            % if stage is not provided, sparse parameter update is performed for all stages.
            obj.t_ocp.set_params_sparse(varargin{:});
        end

        function set_p_global_and_precompute_dependencies(obj, val)
            % usage:
            % ocp.set_p_global_and_precompute_dependencies(val)
            % Sets p_global to val and precomputes all parts of the CasADi graphs of all other functions that only depend on p_global.
            obj.t_ocp.set('p_global', val);
        end


        % function delete(obj)
        %     Use default implementation.
        %     MATLAB destroys the property values after the destruction of the object.
        %     Because `t_ocp` is the only referrence to the `mex_solver` object, MATLAB also destroys the latter.
        % end


    end % methods

    methods (Access = private)
        function generate(obj)

            % generate
            check_dir_and_create(fullfile(pwd, obj.ocp.code_export_directory));
            context = obj.ocp.generate_external_functions();

            obj.ocp.dump_to_json()
            obj.ocp.render_templates()
        end

        function compile_mex_interface_if_needed(obj, solver_creation_opts)

            % check if path contains spaces
            [~,~] = mkdir(solver_creation_opts.output_dir);
            addpath(solver_creation_opts.output_dir);
            if ~isempty(strfind(solver_creation_opts.output_dir, ' '))
                error(strcat('AcadosOcpSolver: Path should not contain spaces, got: ',...
                    solver_creation_opts.output_dir));
            end

            % auto detect whether to compile the interface or not
            if isempty(solver_creation_opts.compile_interface)
                % check if mex interface exists already
                if is_octave()
                    mex_exists = exist( fullfile(solver_creation_opts.output_dir,...
                        '/ocp_get.mex'), 'file');
                else
                    mex_exists = exist( fullfile(solver_creation_opts.output_dir,...
                        ['ocp_get.', mexext]), 'file');
                end
                % check if mex interface is linked against the same external libs as the core
                if mex_exists
                    acados_folder = getenv('ACADOS_INSTALL_DIR');
                    addpath(fullfile(acados_folder, 'external', 'jsonlab'));

                    json_filename = fullfile(acados_folder, 'lib', 'link_libs.json');
                    if ~exist(json_filename, 'file')
                        error('File %s not found.\nPlease compile acados with the latest version, using cmake.', json_filename)
                    end
                    core_links = loadjson(fileread(json_filename));

                    json_filename = fullfile(solver_creation_opts.output_dir, 'link_libs.json');
                    if ~exist(json_filename, 'file')
                        solver_creation_opts.compile_interface = true;
                    else
                        interface_links = loadjson(fileread(json_filename));
                        if isequal(core_links, interface_links)
                            solver_creation_opts.compile_interface = false;
                        else
                            solver_creation_opts.compile_interface = true;
                        end
                    end
                else
                    solver_creation_opts.compile_interface = true;
                end
            end

            if solver_creation_opts.compile_interface
                ocp_compile_interface(solver_creation_opts.output_dir);
                disp('acados MEX interface compiled successfully')
            else
                disp('found compiled acados MEX interface')
            end
        end


        function compile_ocp_shared_lib(self, export_dir)
            return_dir = pwd;
            cd(export_dir);
            if isunix
                %% old code for make
                if ~is_octave()
                    % use Make build system
                    [ status, result ] = system('make ocp_shared_lib');
                    if status
                        cd(return_dir);
                        error('Building templated code as shared library failed.\nGot status %d, result: %s',...
                            status, result);
                    end
                else
                    % use CMake build system, has issues on Linux with MATLAB, see https://github.com/acados/acados/issues/1209
                    [ status, result ] = system('cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_ACADOS_OCP_SOLVER_LIB=ON -S . -B .');
                    if status
                        cd(return_dir);
                        error('Generating buildsystem failed.\nGot status %d, result: %s',...
                            status, result);
                    end
                    [ status, result ] = system('cmake --build . --config Release');
                    if status
                        cd(return_dir);
                        error('Building templated code as shared library failed.\nGot status %d, result: %s',...
                            status, result);
                    end
                end
            else
                % check compiler
                use_msvc = false;
                if ~is_octave()
                    mexOpts = mex.getCompilerConfigurations('C', 'Selected');
                    if contains(mexOpts.ShortName, 'MSVC')
                        use_msvc = true;
                    end
                end
                % compile on Windows platform
                if use_msvc
                    % get env vars for MSVC
                    % msvc_env = fullfile(mexOpts.Location, 'VC\Auxiliary\Build\vcvars64.bat');
                    % assert(isfile(msvc_env), 'Cannot find definition of MSVC env vars.');
                    % detect MSVC version
                    msvc_ver_str = "Visual Studio " + mexOpts.Version(1:2) + " " + mexOpts.Name(22:25);
                    [ status, result ] = system(['cmake -G "' + msvc_ver_str + '" -A x64 -DCMAKE_BUILD_TYPE=Release -DBUILD_ACADOS_OCP_SOLVER_LIB=ON -S . -B .']);
                else
                    [ status, result ] = system('cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DBUILD_ACADOS_OCP_SOLVER_LIB=ON -S . -B .');
                end
                if status
                    cd(return_dir);
                    error('Generating buildsystem failed.\nGot status %d, result: %s',...
                        status, result);
                end
                [ status, result ] = system('cmake --build . --config Release');
                if status
                    cd(return_dir);
                    error('Building templated code as shared library failed.\nGot status %d, result: %s',...
                        status, result);
                end
            end
            cd(return_dir);
        end

    end % private methods

end % class

