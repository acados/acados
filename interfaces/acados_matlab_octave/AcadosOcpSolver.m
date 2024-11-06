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

    properties
        t_ocp % templated solver
        ocp % Matlab class AcadosOcp describing the OCP formulation
        qp_gettable_fields = {'qp_Q', 'qp_R', 'qp_S', 'qp_q', 'qp_r', 'qp_A', 'qp_B', 'qp_b', 'qp_C', 'qp_D', 'qp_lg', 'qp_ug', 'qp_lbx', 'qp_ubx', 'qp_lbu', 'qp_ubu', 'qp_zl', 'qp_zu', 'qp_Zl', 'qp_Zu'}
    end % properties


    methods

        function obj = AcadosOcpSolver(ocp, output_dir)

            if nargin < 2
                output_dir = fullfile(pwd, 'build');
            end

            % detect dimensions & sanity checks
            obj.ocp = ocp;
            obj.ocp.make_consistent()

            % compile mex interface if needed
            obj.compile_mex_interface_if_needed(output_dir)

            % generate
            check_dir_and_create(fullfile(pwd, ocp.code_export_directory));
            context = ocp.generate_external_functions();

            ocp.dump_to_json()
            ocp.render_templates()

            % build
            acados_template_mex.compile_ocp_shared_lib(ocp.code_export_directory)

            % templated MEX
            return_dir = pwd();
            cd(obj.ocp.code_export_directory)

            mex_solver_name = sprintf('%s_mex_solver', obj.ocp.name);
            mex_solver = str2func(mex_solver_name);
            obj.t_ocp = mex_solver();
            addpath(pwd());

            cd(return_dir);
        end


        function solve(obj)
            obj.t_ocp.solve();
        end

        % TODO: remove this! in v.0.5.0!
        function generate_c_code(obj, simulink_opts)
            warning('acados_ocp will be deprecated in the future. Use AcadosOcpSolver instead. For more information on the major acados Matlab interface overhaul, see https://github.com/acados/acados/releases/tag/v0.4.0');
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

        function set(obj, field, value, varargin)
            obj.t_ocp.set(field, value, varargin{:});
        end

        function status = custom_update(obj, data)
            status = obj.t_ocp.custom_update(data);
        end

        function value = get(obj, field, varargin)

            if strcmp('hess_block', field)

                if length(varargin) > 0
                    n = varargin{1};
                    Q = obj.get('qp_Q', n);
                    R = obj.get('qp_R', n);
                    S = obj.get('qp_S', n);

                    value = [R, S; S', Q];
                    return;
                else
                    value = cell(obj.ocp.solver_options.N_horizon, 1);
                    for n=0:obj.ocp.solver_options.N_horizon
                        Q = obj.get('qp_Q', n);
                        R = obj.get('qp_R', n);
                        S = obj.get('qp_S', n);

                        value{n+1} = [R, S; S', Q];

                    end
                    return;
                end
            elseif strcmp('pc_hess_block', field)

                if ~strncmp(obj.ocp.solver_options.qp_solver, 'PARTIAL_CONDENSING', length('PARTIAL_CONDENSING'))
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

        function iterate = get_iterate(obj, iteration)
            if iteration > obj.get('nlp_iter')
                error("iteration needs to be nonnegative and <= nlp_iter.");
            end

            if ~obj.ocp.solver_options.store_iterates
                error("get_iterate: the solver option store_iterates needs to be true in order to get iterates.");
            end

            if strcmp(obj.ocp.solver_options.nlp_solver_type, 'SQP_RTI')
                error("get_iterate: SQP_RTI not supported.");
            end

            N_horizon = obj.ocp.solver_options.N_horizon;

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
            % - min_ev: minimum eigenvalue for each Hessian block.
            % - max_ev: maximum eigenvalue for each Hessian block.
            % - condition_number: condition number for each Hessian block.
            % - condition_number_global: condition number for the full Hessian.

            if length(varargin) > 0
                partially_condensed_qp = varargin{1};
            else
                partially_condensed_qp = false;
            end

            if partially_condensed_qp
                num_blocks = obj.ocp.solver_options.qp_solver_cond_N + 1;
            else
                num_blocks = obj.ocp.dims.N + 1;
            end
            result = struct();
            result.min_ev = zeros(num_blocks, 1);
            result.max_ev = zeros(num_blocks, 1);
            result.condition_number = zeros(num_blocks, 1);
            min_abs_val = 1e12
            max_abs_val = -1e12
            max_ev = -1e12
            min_ev = 1e12

            for n=1:num_blocks
                if partially_condensed_qp
                    hess_block = obj.get('pc_hess_block', n-1);
                else
                    hess_block = obj.get('hess_block', n-1);
                end
                eigvals = eig(hess_block);
                result.min_ev(n) = min(eigvals);
                max_ev = max(max_ev, max(eigvals));
                max_abs_val = max(max_abs_val, max(abs(eigvals)));
                min_ev = min(min_ev, min(eigvals));
                min_abs_val = min(min_abs_val, min(abs(eigvals)));
                result.max_ev(n) = max(eigvals);
                result.condition_number_blockwise(n) = max(eigvals) / min(eigvals);
            end
            result.condition_number_global = max_abs_val / min_abs_val;
            result.max_ev_global = max_ev;
            result.max_abs_ev_global = max_abs_val;
            result.min_ev_global = min_ev;
            result.min_abs_ev_global = min_abs_val;
        end


        function dump_last_qp_to_json(obj, filename)
            qp_data = struct();

            lN = length(num2str(obj.ocp.solver_options.N_horizon+1));
            n_fields = length(obj.qp_gettable_fields);
            for n=1:n_fields

                field = obj.qp_gettable_fields{n};
                for i=0:obj.ocp.solver_options.N_horizon-1
                    s_indx = sprintf(strcat('%0', num2str(lN), 'd'), i);
                    key = strcat(field, '_', s_indx);
                    val = obj.get(field, i);
                    qp_data = setfield(qp_data, key, val);
                end

                if strcmp(field, 'qp_Q') || strcmp(field, 'qp_q') || strcmp(field, 'qp_zl') || strcmp(field, 'qp_zu') || strcmp(field, 'qp_Zl') || strcmp(field, 'qp_Zu')
                    s_indx = sprintf(strcat('%0', num2str(lN), 'd'), obj.ocp.solver_options.N_horizon);
                    key = strcat(field, '_', s_indx);
                    val = obj.get(field, obj.ocp.solver_options.N_horizon);
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

        function compile_mex_interface_if_needed(obj, output_dir)

            % check if path contains spaces
            [~,~] = mkdir(output_dir);
            addpath(output_dir);
            if ~isempty(strfind(output_dir, ' '))
                error(strcat('AcadosOcpSolver: Path should not contain spaces, got: ',...
                    output_dir));
            end

            % auto detect whether to compile the interface or not
            if isempty(obj.ocp.solver_options.compile_interface)
                % check if mex interface exists already
                if is_octave()
                    mex_exists = exist( fullfile(output_dir,...
                        '/ocp_get.mex'), 'file');
                else
                    mex_exists = exist( fullfile(output_dir,...
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

                    json_filename = fullfile(output_dir, 'link_libs.json');
                    if ~exist(json_filename, 'file')
                        obj.ocp.solver_options.compile_interface = true;
                    else
                        interface_links = loadjson(fileread(json_filename));
                        if isequal(core_links, interface_links)
                            obj.ocp.solver_options.compile_interface = false;
                        else
                            obj.ocp.solver_options.compile_interface = true;
                        end
                    end
                else
                    obj.ocp.solver_options.compile_interface = true;
                end
            end

            if obj.ocp.solver_options.compile_interface
                ocp_compile_interface(output_dir);
                disp('acados MEX interface compiled successfully')
            else
                disp('found compiled acados MEX interface')
            end
        end

    end

end % class

