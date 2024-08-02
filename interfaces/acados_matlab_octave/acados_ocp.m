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

classdef acados_ocp < handle

    properties
        % templated solver
        t_ocp
        % matlab objects
        code_gen_dir
        model_struct
        opts_struct
        ocp
        cost_ext_fun_type
        cost_ext_fun_type_e
        cost_ext_fun_type_0
        dyn_ext_fun_type

        qp_gettable_fields = {'qp_Q', 'qp_R', 'qp_S', 'qp_q', 'qp_r', 'qp_A', 'qp_B', 'qp_b', 'qp_C', 'qp_D', 'qp_lg', 'qp_ug', 'qp_lbx', 'qp_ubx', 'qp_lbu', 'qp_ubu'}
    end % properties



    methods


        function obj = acados_ocp(model, opts, simulink_opts)
            obj.model_struct = model.model_struct;
            obj.opts_struct = opts.opts_struct;

            [~,~] = mkdir(obj.opts_struct.output_dir);
            addpath(obj.opts_struct.output_dir);

            % check model consistency
            obj.model_struct = create_consistent_empty_fields(obj.model_struct);

            % detect GNSF structure
            if (strcmp(obj.opts_struct.sim_method, 'irk_gnsf'))
                if (strcmp(obj.opts_struct.gnsf_detect_struct, 'true'))
                    obj.model_struct = detect_gnsf_structure(obj.model_struct);
                    generate_get_gnsf_structure(obj.model_struct, obj.opts_struct);
                else
                    obj.model_struct = get_gnsf_structure(obj.model_struct);
                end
            end

            % store ext_fun_type
            obj.cost_ext_fun_type = obj.model_struct.cost_ext_fun_type;
            obj.cost_ext_fun_type_e = obj.model_struct.cost_ext_fun_type_e;
            obj.cost_ext_fun_type_0 = obj.model_struct.cost_ext_fun_type_0;
            obj.dyn_ext_fun_type = obj.model_struct.dyn_ext_fun_type;

            % detect cost type
            if (strcmp(obj.model_struct.cost_type, 'auto'))
                obj.model_struct = detect_cost_type(obj.model_struct, 'path');
            end
            if (strcmp(obj.model_struct.cost_type_0, 'auto'))
                obj.model_struct = detect_cost_type(obj.model_struct, 'initial');
            elseif isempty(obj.model_struct.cost_type_0)
                % copy entries from path cost
                obj.model_struct.cost_type_0 = obj.model_struct.cost_type;
                if (strcmp(obj.model_struct.cost_type, 'linear_ls'))
                    obj.model_struct.cost_Vx_0 = obj.model_struct.cost_Vx;
                    obj.model_struct.cost_Vu_0 = obj.model_struct.cost_Vu;
                    if isfield(obj.model_struct, 'cost_Vz')
                        obj.model_struct.cost_Vz_0 = obj.model_struct.cost_Vz;
                    end
                elseif (strcmp(obj.model_struct.cost_type, 'nonlinear_ls'))
                    obj.model_struct.cost_expr_y_0 = obj.model_struct.cost_expr_y;
                elseif (strcmp(obj.model_struct.cost_type, 'ext_cost'))
                    obj.model_struct.cost_ext_fun_type_0 = obj.model_struct.cost_ext_fun_type;
                    if strcmp(obj.model_struct.cost_ext_fun_type_0, 'casadi')
                        obj.model_struct.cost_expr_ext_cost_0 = obj.model_struct.cost_expr_ext_cost;
                        if isfield(obj.model_struct, 'cost_expr_ext_cost_custom_hess')
                            obj.model_struct.cost_expr_ext_cost_custom_hess_0 = obj.model_struct.cost_expr_ext_cost_custom_hess;
                        end
                    else % generic
                        obj.model_struct.cost_source_ext_cost_0 = obj.model_struct.cost_source_ext_cost;
                        obj.model_struct.cost_function_ext_cost_0 = obj.model_struct.cost_function_ext_cost;
                    end
                end
                if (strcmp(obj.model_struct.cost_type, 'linear_ls')) || (strcmp(obj.model_struct.cost_type, 'nonlinear_ls'))
                    obj.model_struct.cost_W_0 = obj.model_struct.cost_W;
                    if isfield(obj.model_struct,'cost_y_ref')
                        obj.model_struct.cost_y_ref_0 = obj.model_struct.cost_y_ref;
                    end
                end
            end
            if (strcmp(obj.model_struct.cost_type_e, 'auto'))
                obj.model_struct = detect_cost_type(obj.model_struct, 'terminal');
            end

            % detect constraint structure
            if (strcmp(obj.model_struct.constr_type, 'auto'))
                obj.model_struct = detect_constr(obj.model_struct, 0);
            end
            if (strcmp(obj.model_struct.constr_type_0, 'auto'))
                obj.model_struct = detect_constr(obj.model_struct, 0, true);
            end
            if (strcmp(obj.model_struct.constr_type_e, 'auto'))
                obj.model_struct = detect_constr(obj.model_struct, 1);
            end

            % detect dimensions & sanity checks
            [obj.model_struct, obj.opts_struct] = detect_dims_ocp(obj.model_struct, obj.opts_struct);

            % check if path contains spaces
            if ~isempty(strfind(obj.opts_struct.output_dir, ' '))
                error(strcat('acados_ocp: Path should not contain spaces, got: ',...
                    obj.opts_struct.output_dir));
            end

            % compile mex interface (without model dependency)
            if ( strcmp(obj.opts_struct.compile_interface, 'true') )
                compile_interface = true;
            elseif ( strcmp(obj.opts_struct.compile_interface, 'false') )
                compile_interface = false;
            elseif ( strcmp(obj.opts_struct.compile_interface, 'auto') )
                % check if mex interface exists already
                if is_octave()
                    mex_exists = exist( fullfile(obj.opts_struct.output_dir,...
                        '/ocp_get.mex'), 'file');
                else
                    mex_exists = exist( fullfile(obj.opts_struct.output_dir,...
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

                    json_filename = fullfile(obj.opts_struct.output_dir, 'link_libs.json');
                    if ~exist(json_filename, 'file')
                        compile_interface = true;
                    else
                        interface_links = loadjson(fileread(json_filename));
                        if isequal(core_links, interface_links)
                            compile_interface = false;
                        else
                            compile_interface = true;
                        end
                    end
                else
                    compile_interface = true;
                end
            else
                error('acados_ocp: field compile_interface is %, supported values are: true, false, auto', ...
                        obj.opts_struct.compile_interface);
            end

            if ( compile_interface )
                ocp_compile_interface(obj.opts_struct);
                disp('acados MEX interface compiled successfully')
            else
                disp('found compiled acados MEX interface')
            end

            % check for unsupported options:
            if strcmp(obj.opts_struct.qp_solver, "partial_condensing_osqp") || strcmp(obj.opts_struct.qp_solver, "partial_condensing_hpmpc") || strcmp(obj.opts_struct.qp_solver, "partial_condensing_qpdunes") || ...
                strcmp(obj.opts_struct.qp_solver, "partial_condensing_ooqp")
                if obj.model_struct.dim_ns > 0 || obj.model_struct.dim_ns_e > 0
                    error(['selected QP solver ', obj.opts_struct.qp_solver, ' does not support soft constraints (yet).'])
                end
            end

            % generate templated solver
            if nargin < 3
                simulink_opts = get_acados_simulink_opts();
            end
            obj.ocp = setup_ocp(obj, simulink_opts);
            ocp_generate_c_code(obj.ocp);

            % templated MEX
            return_dir = pwd();
            obj.code_gen_dir = obj.ocp.code_export_directory;
            cd(obj.code_gen_dir)

            mex_solver_name = sprintf('%s_mex_solver', obj.model_struct.name);
            mex_solver = str2func(mex_solver_name);
            obj.t_ocp = mex_solver();
            addpath(pwd());

            cd(return_dir);

        end


        function solve(obj)
            obj.t_ocp.solve();
        end

        % TODO: remove this? does not seem to do anything
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
            if strcmp('qp_Q', field) || strcmp('qp_R', field) || strcmp('qp_solver_cond_H', field)
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

        function result = qp_diagnostics(obj, varargin)
            % Compute some diagnostic values for the last QP.
            % min_ev: minimum eigenvalue for each Hessian block.
            % max_ev: maximum eigenvalue for each Hessian block.
            % condition_number: condition number for each Hessian block.

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

            for n=1:num_blocks
                if partially_condensed_qp
                    disp(['pc_hess_block', num2str(n-1)]);
                    hess_block = obj.get('pc_hess_block', n-1);
                else
                    hess_block = obj.get('hess_block', n-1);
                end

                eigvals = eig(hess_block);
                result.min_ev(n) = min(eigvals);
                result.max_ev(n) = max(eigvals);
                result.condition_number_blockwise(n) = max(eigvals) / min(eigvals);
            end
            result.condition_number_global = max(result.max_ev) / min(result.min_ev);
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

end % class

