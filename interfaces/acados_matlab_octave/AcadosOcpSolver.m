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
        solver_creation_opts
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
            % - check_reuse_possible: boolean, default true.
            %        if true and generate is false:
            %        check if code reuse is possible by comparing OCP formulations,
            %        options and acados version. If not identical, code generation and build are forced.
            % - compile_mex_wrapper: boolean, if true, the mex wrapper is compiled
            % - compile_interface: can be [], true or false. If [], the interface is compiled if it does not exist.
            % - output_dir: path to the directory where the MEX interface is compiled
            obj.ocp = ocp;

            % optional arguments
            % solver creation options
            default_solver_creation_opts = struct('json_file', '', ...
                    'build', true, ...
                    'generate', true, ...
                    'check_reuse_possible', true, ...
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
            obj.solver_creation_opts = solver_creation_opts;

            if isempty(ocp) && isempty(solver_creation_opts.json_file)
                error('AcadosOcpSolver: provide either an OCP object or a json file');
            end

            if isempty(ocp)
                json_file = obj.solver_creation_opts.json_file;
                if obj.solver_creation_opts.generate
                    disp('AcadosOcpSolver: OCP not provided, cannot generate code, setting generate to false');
                    obj.solver_creation_opts.generate = false;
                end
                obj.solver_creation_opts.check_reuse_possible = false;
            else
                % formulation provided
                if ~isempty(ocp.solver_options.compile_interface) && ~isempty(obj.solver_creation_opts.compile_interface)
                    error('AcadosOcpSolver: provide either compile_interface in OCP object or obj.solver_creation_opts');
                end
                if ~isempty(ocp.solver_options.compile_interface)
                    obj.solver_creation_opts.compile_interface = ocp.solver_options.compile_interface;
                end
                if ~isempty(obj.solver_creation_opts.json_file)
                    ocp.code_gen_opts.json_file = obj.solver_creation_opts.json_file;
                end
                % make consistent
                ocp.make_consistent();

                json_file = ocp.code_gen_opts.json_file;
            end

            %% compile mex interface if needed
            obj.compile_mex_interface_if_needed();

            %% generate
            if ~obj.solver_creation_opts.generate && obj.solver_creation_opts.check_reuse_possible
                % check if code reuse can be done
                reuse_possible = obj.is_code_reuse_possible(json_file, 1);
                if ~reuse_possible
                    disp('AcadosOcpSolver: code reuse not possible, forcing code generation and build...');
                    obj.solver_creation_opts.generate = true;
                    obj.solver_creation_opts.build = true;
                else
                    disp('AcadosOcpSolver: attempting code reuse...')
                end
            end

            if obj.solver_creation_opts.generate
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
            code_export_directory = acados_ocp_struct.code_gen_opts.code_export_directory;

            %% compile problem specific shared library
            if obj.solver_creation_opts.build
                obj.compile_ocp_shared_lib(code_export_directory);
            end

            %% create solver
            return_dir = pwd();
            cd(code_export_directory)

            mex_solver_name = sprintf('%s_mex_solver', obj.name);
            mex_solver = str2func(mex_solver_name);
            obj.t_ocp = mex_solver(obj.solver_creation_opts);
            addpath(pwd());

            cd(return_dir);
        end

        function code_reuse_possible = is_code_reuse_possible(obj, json_file, verbose)
            code_reuse_possible = 1;
            if ~exist(obj.ocp.code_gen_opts.code_export_directory, 'dir')
                code_reuse_possible = 0;
                if verbose
                    disp('code reuse not possible: code export directory does not exist');
                end
                return;
            end
            if ~exist(json_file, 'file')
                code_reuse_possible = 0;
                if verbose
                    disp('code reuse not possible: json file does not exist');
                end
                return;
            end
            try
                ocp_struct_restore = loadjson(fileread(json_file), 'SimplifyCell', 0);
            catch
                code_reuse_possible = 0;
                if verbose
                    disp('code reuse not possible: error loading json file');
                end
                return;
            end

            try
                old_hash = ocp_struct_restore.hash;
            catch
                code_reuse_possible = 0;
                if verbose
                    disp('code reuse not possible: no hash in json file');
                end
                return;
            end

            % disp(['old hash', ocp_struct_restore.hash])
            % old_ocp = AcadosOcp.from_struct(ocp_struct_restore);

            % create hash for current ocp
            try
                obj.ocp.make_consistent()
                ocp_struct = orderfields(obj.ocp.to_struct());
                new_hash = hash_struct(ocp_struct);
            catch
                code_reuse_possible = 0;
                if verbose
                    disp('code reuse not possible: error creating hash for current ocp');
                end
                return;
            end

            if strcmp(old_hash, new_hash) ~= 1
                code_reuse_possible = 0;
                if verbose
                    disp('code reuse not possible: hash mismatch');
                end
                return;
            end
        end

        function solve(obj)
            obj.t_ocp.solve();
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
            if strcmp(obj.solver_options.qpscaling_scale_constraints, 'NO_CONSTRAINT_SCALING')
                error("QP scaling constraints are not enabled.");
            end
            % returns the qp scaling constraints for the given stage
            value = obj.t_ocp.get('qpscaling_constr', stage);
        end

        function value = get_qp_scaling_objective(obj)
            % returns the qp scaling objective
            if strcmp(obj.solver_options.qpscaling_scale_objective, 'NO_OBJECTIVE_SCALING')
                error("QP scaling objective is not enabled.");
            end
            value = obj.t_ocp.get('qpscaling_obj');
        end

        function value = get(obj, field, varargin)
            if strcmp('qp_iter_all', field)
                full_stats = obj.t_ocp.get('stat');
                if strcmp(obj.solver_options.nlp_solver_type, 'SQP')
                    value = full_stats(:, 7);
                    return;
                elseif strcmp(obj.solver_options.nlp_solver_type, 'SQP_RTI')
                    value = full_stats(:, 3);
                    return;
                else
                    error("qp_iter is not available for nlp_solver_type %s.", obj.solver_options.nlp_solver_type);
                end
            end

            if strcmp('res_all', field)
                if ~strcmp(obj.solver_options.nlp_solver_type, 'SQP')
                    error("res_all is only available for nlp_solver_type SQP.");
                end
                full_stats = obj.t_ocp.get('stat');
                value = full_stats(:, 2:5);
                return;
            end

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
            filename = '';
            overwrite = false;

            if nargin>=2
                filename = varargin{1};
                if ~isa(filename, 'char')
                    error('filename must be a char vector, use '' ''');
                end
            end

            if nargin==3
                overwrite = varargin{2};
            end

            if nargin > 3
                disp('acados_ocp.get: wrong number of input arguments (1 or 2 allowed)');
            end

            if strcmp(filename,'')
                filename = [obj.name '_iterate.json'];
            end
            if ~overwrite
                % append timestamp
                if exist(filename, 'file')
                    filename = filename(1:end-5);
                    filename = [filename '_' datestr(now,'yyyy-mm-dd-HH:MM:SS') '.json'];
                end
            end
            filename = fullfile(pwd, filename);

            % get iterate:
            solution = struct();
            for i=0:obj.N_horizon
                solution.(['x_' num2str(i)]) = obj.get('x', i);
                solution.(['lam_' num2str(i)]) = obj.get('lam', i);
                solution.(['sl_' num2str(i)]) = obj.get('sl', i);
                solution.(['su_' num2str(i)]) = obj.get('su', i);
            end
            for i=0:obj.N_horizon-1
                solution.(['z_' num2str(i)]) = obj.get('z', i);
                solution.(['u_' num2str(i)]) = obj.get('u', i);
                solution.(['pi_' num2str(i)]) = obj.get('pi', i);
            end

            acados_folder = getenv('ACADOS_INSTALL_DIR');
            addpath(fullfile(acados_folder, 'external', 'jsonlab'));
            savejson('', solution, filename);

            json_string = savejson('', solution, 'ForceRootName', 0);

            fid = fopen(filename, 'w');
            if fid == -1, error('store_iterate: Cannot create JSON file'); end
            fwrite(fid, json_string, 'char');
            fclose(fid);

            disp(['stored current iterate in ' filename]);
        end


        function [] = load_iterate(obj, filename)
            %%%  Loads the iterate stored in json file with filename into the ocp solver.
            acados_folder = getenv('ACADOS_INSTALL_DIR');
            addpath(fullfile(acados_folder, 'external', 'jsonlab'));
            filename = fullfile(pwd, filename);

            if ~exist(filename, 'file')
                error(['load_iterate: failed, file does not exist: ' filename])
            end

            solution = loadjson(filename);
            keys = fieldnames(solution);

            for k = 1:numel(keys)
                key = keys{k};
                key_parts = strsplit(key, '_');
                field = key_parts{1};
                stage = key_parts{2};

                val = solution.(key);

                % check if array is empty (can happen for z)
                if numel(val) > 0
                    obj.set(field, val, str2num(stage))
                end
            end
        end

        function iterate = store_iterate_to_obj(obj)
            % Returns the current iterate of the OCP solver as an AcadosOcpIterate.

            warning('Deprecation warning: store_iterate_to_obj() is deprecated, use get_iterate() instead.');
            iterate = obj.get_iterate();
        end

        function [] = load_iterate_from_obj(obj, iterate)
            %%%  Loads the iterate from an AcadosOcpIterate object.
            %%% param1: iterate: AcadosOcpIterate object containing the iterate to load
            warning('Deprecation warning: load_iterate_from_obj() is deprecated, use set_iterate() instead.');
            obj.set_iterate(iterate);
        end

        function [] = set_iterate(obj, iterate)
            %%%  Loads the iterate from an AcadosOcpIterate object.
            %%% param1: iterate: AcadosOcpIterate object containing the iterate to load

            if ~isa(iterate, 'AcadosOcpIterate')
                error('set_iterate: iterate needs to be of type AcadosOcpIterate');
            end

            for i = 1:obj.solver_options.N_horizon + 1
                obj.t_ocp.set('x', iterate.x{i, 1}, i-1);
                obj.t_ocp.set('sl', iterate.sl{i, 1}, i-1);
                obj.t_ocp.set('su', iterate.su{i, 1}, i-1);
                obj.t_ocp.set('lam', iterate.lam{i, 1}, i-1);
            end
            for i = 1:obj.solver_options.N_horizon
                obj.t_ocp.set('u', iterate.u{i, 1}, i-1);
                obj.t_ocp.set('pi', iterate.pi{i, 1}, i-1);
                if ~isempty(iterate.z{i, 1})
                    obj.t_ocp.set('z', iterate.z{i, 1}, i-1);
                end
            end
        end
        function iterate = get_iterate(obj, iteration)
            nlp_iter = obj.get('nlp_iter');

            get_last_iterate = nargin == 1 || iteration == -1 || iteration == nlp_iter;

            if ~get_last_iterate && (iteration > nlp_iter || iteration < 0)
                error("iteration needs to be nonnegative and <= nlp_iter.");
            end

            if ~get_last_iterate && ~obj.solver_options.store_iterates
                error("get_iterate: the solver option store_iterates needs to be true in order to get intermediate iterates.");
            end

            if ~get_last_iterate && strcmp(obj.solver_options.nlp_solver_type, 'SQP_RTI')
                error("get_iterate: SQP_RTI not supported.");
            end

            fields = {'x','u','z','sl','su','pi','lam'};
            d = struct();
            for fi = 1:length(fields)
                field = fields{fi};
                traj = {};
                for n = 0:obj.solver_options.N_horizon
                    if n < obj.solver_options.N_horizon || ~ismember(field, {'u','pi','z'})

                        if get_last_iterate
                            val = obj.get(field, n);
                        else
                            val = obj.get(field, n, iteration);
                        end
                        traj{end+1,1} = val;
                    end
                end
                d.(field) = traj;
            end

            iterate = AcadosOcpIterate( ...
                d.x, d.u, d.z, ...
                d.sl, d.su, d.pi, d.lam );
        end

        function iterates = get_iterates(obj)
            nlp_iter = obj.get('nlp_iter');
            iterates_cell = cell(nlp_iter+1, 1);

            for n=1:(nlp_iter+1)
                iterates_cell{n} = obj.get_iterate(n-1);
            end

            iterates = AcadosOcpIterates(iterates_cell);
        end

        function qp_data = get_last_qp(obj)
            %%% Returns the latest QP data as a struct
            qp_data = struct();

            % Define QP field groups similar to Python implementation
            qp_dynamics_fields = {'A', 'B', 'b'};
            qp_cost_fields = {'Q', 'R', 'S', 'q', 'r', 'zl', 'zu', 'Zl', 'Zu'};
            qp_constraint_fields = {'C', ...
                                'D', ...
                                'lg', ...
                                'ug', ...
                                'lbx', ...
                                'ubx', ...
                                'lbu', ...
                                'ubu', ...
                                'lls', ...
                                'lus', ...
                                'lg_mask', ...
                                'ug_mask', ...
                                'lbx_mask', ...
                                'ubx_mask', ...
                                'lbu_mask', ...
                                'ubu_mask', ...
                                'lls_mask', ...
                                'lus_mask', ...
                                'idxs_rev', ...
                                'idxb', ...
                                'idxe'};

            % Format for zero-padding stage numbers
            lN = length(num2str(obj.N_horizon));

            % Get cost and constraint fields (for stages 0 to N)
            all_other_fields = [qp_cost_fields, qp_constraint_fields, qp_dynamics_fields];
            for i = 1:length(all_other_fields)
                field = all_other_fields{i};
                k_last = obj.N_horizon;
                if ismember(field, qp_dynamics_fields)
                    k_last = obj.N_horizon - 1; % dynamics only up to N-1
                end
                for stage = 0:k_last
                    field_name = sprintf('%s_%0*d', field, lN, stage);
                    value = obj.get(strcat('qp_', field), stage);
                    if ~isempty(value)
                        qp_data.(field_name) = value;
                    end
                end
            end
        end

        function dump_last_qp_to_json(obj, varargin)
            %%% Dumps the latest QP data into a json file
            %%% param1: filename: if not set, use model_name + '_QP.json'
            %%% param2: overwrite: if false and filename exists add timestamp to filename
            %%% param3: backend: 'Matlab' (default) or 'C'.
            filename = '';
            overwrite = false;

            if nargin>=2
                filename = varargin{1};
                if nargin>=3
                    overwrite = varargin{2};
                end
                backend = 'Matlab';
                if nargin >= 4
                    backend = varargin{3};
                end
            end

            if nargin > 4
                disp('dump_last_qp_to_json: wrong number of input arguments (1 or 2 allowed)');
            end

            if strcmp(filename,'')
                filename = [obj.name '_QP.json'];
            end
            if ~overwrite
                if exist(filename, 'file')
                    % append timestamp
                    [~, name, ~] = fileparts(filename);
                    timestamp = datestr(now,'yyyy-mm-dd-HH-MM-SS.FFF');
                    filename = [name '_' timestamp '.json'];
                end
            end
            filename = fullfile(pwd, filename);

            if strcmp(backend, 'Matlab')
                % get QP data:
                qp_data = obj.get_last_qp();

                acados_folder = getenv('ACADOS_INSTALL_DIR');
                addpath(fullfile(acados_folder, 'external', 'jsonlab'));

                json_string = savejson('', qp_data, 'ForceRootName', 0);

                fid = fopen(filename, 'w');
                if fid == -1, error('dump_last_qp_to_json: Cannot create JSON file'); end
                fwrite(fid, json_string, 'char');
                fclose(fid);

                disp(['stored qp from solver memory in ' filename]);
            elseif strcmp(backend, 'C')
                obj.t_ocp.dump_last_qp_to_json(filename);
                disp(['stored qp with C backend from solver memory in ' filename]);
             else
                error('dump_last_qp_to_json: backend not recognized, use ''Matlab'' or ''C''.');
            end
        end

        function print(obj, varargin)
            obj.t_ocp.print(varargin{:});
        end

        function print_statistics(obj)
            obj.t_ocp.print();
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
            max_abs_eig_val = -inf;
            max_ev = -inf;
            min_ev = inf;
            max_abs_hess_val = -inf;

            for n=1:num_blocks
                if partially_condensed_qp
                    hess_block = obj.get('pc_hess_block', n-1);
                else
                    hess_block = obj.get('hess_block', n-1);
                end
                eigvals = eig(hess_block);
                max_ev = max(max_ev, max(eigvals));
                max_abs_eig_val = max(max_abs_eig_val, max(abs(eigvals)));
                min_ev = min(min_ev, min(eigvals));
                min_abs_val = min(min_abs_val, min(abs(eigvals)));
                max_abs_hess_val = max(max_abs_hess_val, max(max(abs(hess_block))));
                result.min_eigv_stage(n) = min(eigvals);
                result.max_eigv_stage(n) = max(eigvals);
                result.condition_number_stage(n) = max(eigvals) / min(eigvals);
            end
            % Zl, Zu don't change in partial condensing.
            for n=1:obj.N_horizon+1
                max_abs_hess_val = max(max_abs_hess_val, max(abs(obj.get('qp_Zl', n-1))));
                max_abs_hess_val = max(max_abs_hess_val, max(abs(obj.get('qp_Zu', n-1))));
            end
            result.condition_number_global = max_abs_eig_val / min_abs_val;
            result.max_eigv_global = max_ev;
            result.max_abs_eigv_global = max_abs_eig_val;
            result.min_eigv_global = min_ev;
            result.min_abs_eigv_global = min_abs_val;
            result.max_abs_hess_val = max_abs_hess_val;
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
            check_dir_and_create(obj.ocp.code_gen_opts.code_export_directory);
            context = obj.ocp.generate_external_functions();

            obj.ocp.dump_to_json()
            obj.ocp.render_templates()
        end

        function compile_mex_interface_if_needed(obj)

            % check if path contains spaces
            [~,~] = mkdir(obj.solver_creation_opts.output_dir);
            addpath(obj.solver_creation_opts.output_dir);
            if ~isempty(strfind(obj.solver_creation_opts.output_dir, ' '))
                error(strcat('AcadosOcpSolver: Path should not contain spaces, got: ',...
                    obj.solver_creation_opts.output_dir));
            end

            % auto detect whether to compile the interface or not
            if isempty(obj.solver_creation_opts.compile_interface)
                % check if mex interface exists already
                if is_octave()
                    mex_exists = exist( fullfile(obj.solver_creation_opts.output_dir,...
                        '/ocp_get.mex'), 'file');
                else
                    mex_exists = exist( fullfile(obj.solver_creation_opts.output_dir,...
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

                    json_filename = fullfile(obj.solver_creation_opts.output_dir, 'link_libs.json');
                    if ~exist(json_filename, 'file')
                        obj.solver_creation_opts.compile_interface = true;
                    else
                        interface_links = loadjson(fileread(json_filename));
                        if isequal(core_links, interface_links)
                            obj.solver_creation_opts.compile_interface = false;
                        else
                            obj.solver_creation_opts.compile_interface = true;
                        end
                    end
                else
                    obj.solver_creation_opts.compile_interface = true;
                end
            end

            if obj.solver_creation_opts.compile_interface
                ocp_compile_interface(obj.solver_creation_opts.output_dir);
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
            disp('problem specific shared library compiled successfully');
        end

    end % private methods

end % class

