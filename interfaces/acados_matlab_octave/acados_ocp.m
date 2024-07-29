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
        ocp
        cost_ext_fun_type
        cost_ext_fun_type_e
        cost_ext_fun_type_0
        dyn_ext_fun_type
    end % properties



    methods

        function obj = acados_ocp(model, opts, simulink_opts)

            if nargin < 3
                simulink_opts = get_acados_simulink_opts();
            end
            % TODO where do we get these from
            output_dir = opts.opts_struct.output_dir;
            gnsf_transcription_opts = struct();

            [~,~] = mkdir(output_dir);
            addpath(output_dir);

            obj.ocp = setup_ocp(obj, model.model_struct, opts.opts_struct, simulink_opts);
            % detect dimensions & sanity checks
            obj.ocp.model.make_consistent(obj.ocp.dims);
            detect_dims_ocp(obj.ocp);
            obj.ocp.make_consistent()

            % detect GNSF structure
            if strcmp(obj.ocp.solver_options.integrator_type, 'GNSF')
                if obj.ocp.dims.gnsf_nx1 + dims.gnsf_nx2 ~= ocp.dims.nx
                    detect_gnsf_structure(obj.ocp.model, obj.ocp.dims, gnsf_transcription_opts);
                else
                    warning('No GNSF model detected, assuming all required fields are set.')
                end
            end

            % detect cost type
            stage_types = {'initial', 'path', 'terminal'};
            cost_types = {obj.ocp.cost.cost_type_0, obj.ocp.cost.cost_type, obj.ocp.cost.cost_type_e};

            for n=1:3
                if strcmp(cost_types{n}, 'AUTO')
                    detect_cost_type(obj.ocp.model, obj.ocp.cost, obj.ocp.dims, stage_types{n});
                end
            end

            % if initial is empty, copy path cost
            if isempty(cost_types{1})
                warning("cost_type_0 not set, using path cost");
                obj.ocp.cost.cost_type_0 = obj.ocp.cost.cost_type;
                if (strcmp(obj.ocp.cost.cost_type, 'LINEAR_LS'))
                    obj.ocp.cost.Vx_0 = obj.ocp.cost.Vx;
                    obj.ocp.cost.Vu_0 = obj.ocp.cost.Vu;
                    obj.ocp.cost.Vz_0 = obj.ocp.cost.Vz;
                elseif (strcmp(obj.ocp.cost.cost_type, 'NONLINEAR_LS'))
                    obj.ocp.model.cost_y_expr_0 = obj.ocp.model.cost_y_expr;
                elseif (strcmp(obj.ocp.cost.cost_type, 'EXTERNAL'))
                    obj.ocp.cost.cost_ext_fun_type_0 = obj.ocp.cost.cost_ext_fun_type;
                    if strcmp(obj.ocp.cost.cost_ext_fun_type_0, 'casadi')
                        obj.ocp.model.cost_expr_ext_cost_0 = obj.ocp.model.cost_expr_ext_cost;
                        obj.ocp.model.cost_expr_ext_cost_custom_hess_0 = obj.ocp.model.cost_expr_ext_cost_custom_hess;
                    else % generic
                        obj.ocp.cost.cost_source_ext_cost_0 = obj.ocp.cost.cost_source_ext_cost;
                        obj.ocp.cost.cost_function_ext_cost_0 = obj.ocp.cost.cost_function_ext_cost;
                    end
                end
                if (strcmp(obj.ocp.cost.cost_type, 'LINEAR_LS')) || (strcmp(obj.ocp.cost.cost_type, 'NONLINEAR_LS'))
                    obj.ocp.cost.W_0 = obj.ocp.cost.W;
                    obj.ocp.cost.yref_0 = obj.ocp.cost.yref;
                    obj.ocp.dims.ny_0 = obj.ocp.dims.ny;
                end
            end

            % detect constraint structure
            constraint_types = {obj.ocp.constraints.constr_type_0, obj.ocp.constraints.constr_type, obj.ocp.constraints.constr_type_e};
            for n=1:3
                if strcmp(constraint_types{n}, 'AUTO')
                    detect_constr(obj.ocp.model, obj.ocp.constraints, obj.ocp.dims, stage_types{n});
                end
            end


            % check if path contains spaces
            if ~isempty(strfind(output_dir, ' '))
                error(strcat('acados_ocp: Path should not contain spaces, got: ',...
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


            % check for unsupported options:
            if strcmp(obj.ocp.solver_options.nlp_solver_type, "PARTIAL_CONDENSING_OSQP") || ...
                strcmp(obj.ocp.solver_options.nlp_solver_type, "PARTIAL_CONDENSING_HPMPC") || ...
                strcmp(obj.ocp.solver_options.nlp_solver_type, "PARTIAL_CONDENSING_QPDUNES") || ...
                strcmp(obj.ocp.solver_options.nlp_solver_type, "PARTIAL_CONDENSING_OOQP")
                if obj.ocp.dims.ns > 0 || obj.ocp.dims.ns_e > 0
                    error(['selected QP solver ', obj.ocp.solver_options.nlp_solver_type, ' does not support soft constraints (yet).'])
                end
            end

            % generate templated solver
            ocp_generate_c_code(obj.ocp);

            % templated MEX
            return_dir = pwd();
            obj.code_gen_dir = obj.ocp.code_export_directory;
            cd(obj.code_gen_dir)

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
            value = obj.t_ocp.get(field, varargin{:});
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

        % function delete(obj)
        %     Use default implementation.
        %     MATLAB destroys the property values after the destruction of the object.
        %     Because `t_ocp` is the only referrence to the `mex_solver` object, MATLAB also destroys the latter.
        % end


    end % methods

end % class

