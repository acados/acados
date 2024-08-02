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
        qp_gettable_fields = {'qp_Q', 'qp_R', 'qp_S', 'qp_q', 'qp_r', 'qp_A', 'qp_B', 'qp_b', 'qp_C', 'qp_D', 'qp_lg', 'qp_ug', 'qp_lbx', 'qp_ubx', 'qp_lbu', 'qp_ubu'}
    end % properties


    methods

        function obj = AcadosOcpSolver(ocp, output_dir, simulink_opts)

            if nargin < 3
                simulink_opts = get_acados_simulink_opts();
            end

            if nargin < 2
                output_dir = fullfile(pwd, 'build');
            end

            % detect dimensions & sanity checks
            obj.ocp = ocp;
            obj.ocp.make_consistent()

            % detect GNSF structure
            if strcmp(obj.ocp.solver_options.integrator_type, 'GNSF')
                if obj.ocp.dims.gnsf_nx1 + obj.ocp.dims.gnsf_nx2 ~= obj.ocp.dims.nx
                    % TODO: properly interface those.
                    gnsf_transcription_opts = struct();
                    detect_gnsf_structure(obj.ocp.model, obj.ocp.dims, gnsf_transcription_opts);
                else
                    warning('No GNSF model detected, assuming all required fields are set.')
                end
            end

            obj.compile_mex_interface_if_needed(output_dir)

            % generate templated solver
            ocp_generate_c_code(obj.ocp);

            % templated MEX
            return_dir = pwd();
            cd(obj.ocp.code_export_directory)

            mex_solver_name = sprintf('%s_mex_solver', obj.ocp.model.name);
            mex_solver = str2func(mex_solver_name);
            obj.t_ocp = mex_solver();
            addpath(pwd());

            cd(return_dir);
        end


        function solve(obj)
            obj.t_ocp.solve();
        end

        % TODO: remove this? does not seem to do anything
        % remove from all examples first!
        function generate_c_code(obj, simulink_opts)
            if nargin < 2
                warning("Code is generated with the default simulink options via the constructor of acados_ocp.")
            else
                error("If you want to provide simulink options, put it in the constructor of acados_ocp.")
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
                    value = cell(obj.ocp.dims.N, 1);
                    for n=0:obj.ocp.dims.N
                        Q = obj.get('qp_Q', n);
                        R = obj.get('qp_R', n);
                        S = obj.get('qp_S', n);

                        value{n+1} = [R, S; S', Q];

                    end
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
            obj.t_ocp.store_iterate(varargin{:});
        end


        function [] = load_iterate(obj, filename)
            obj.t_ocp.load_iterate(filename);
        end


        function print(obj, varargin)
            obj.t_ocp.print(varargin{:});
        end


        function reset(obj)
            obj.t_ocp.reset();
        end


        function result = qp_diagnostics(obj)
            % Compute some diagnostic values for the last QP.
            % min_ev: minimum eigenvalue for each Hessian block.
            % max_ev: maximum eigenvalue for each Hessian block.
            % condition_number: condition number for each Hessian block.

            result = struct();
            result.min_ev = zeros(obj.ocp.dims.N+1, 1);
            result.max_ev = zeros(obj.ocp.dims.N+1, 1);
            result.condition_number = zeros(obj.ocp.dims.N+1, 1);

            for n=0:obj.ocp.dims.N
                if n < obj.ocp.dims.N
                    Q = obj.get('qp_Q', n);
                    R = obj.get('qp_R', n);
                    S = obj.get('qp_S', n);
                    hess_block = [R, S; S', Q];
                else
                    hess_block = Q;
                end

                eigvals = eig(hess_block);
                result.min_ev(n+1) = min(eigvals);
                result.max_ev(n+1) = max(eigvals);
                result.condition_number(n+1) = max(eigvals) / min(eigvals);
            end
        end


        function dump_last_qp_to_json(obj, filename)
            qp_data = struct();

            lN = length(num2str(obj.ocp.dims.N+1));
            n_fields = length(obj.qp_gettable_fields);
            for n=1:n_fields

                field = obj.qp_gettable_fields{n};
                for i=0:obj.ocp.dims.N-1
                    s_indx = sprintf(strcat('%0', num2str(lN), 'd'), i);
                    key = strcat(field, '_', s_indx);
                    val = obj.get(field, i);
                    qp_data = setfield(qp_data, key, val);
                end

                if strcmp(field, 'qp_Q') || strcmp(field, 'qp_q')
                    s_indx = sprintf(strcat('%0', num2str(lN), 'd'), obj.ocp.dims.N);
                    key = strcat(field, '_', s_indx);
                    val = obj.get(field, obj.ocp.dims.N);
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

