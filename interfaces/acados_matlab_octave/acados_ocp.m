%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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
        C_ocp
        C_ocp_ext_fun
        model_struct
        opts_struct
        acados_ocp_nlp_json
        cost_ext_fun_type
        cost_ext_fun_type_e
        cost_ext_fun_type_0
        dyn_ext_fun_type
    end % properties



    methods


        function obj = acados_ocp(model, opts)
            obj.model_struct = model.model_struct;
            obj.opts_struct = opts.opts_struct;

            [~,~] = mkdir(obj.opts_struct.output_dir);
            addpath(obj.opts_struct.output_dir);

            % check model consistency
            obj.model_struct = create_consistent_empty_fields(obj.model_struct, obj.opts_struct);

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
                        '/ocp_create.mex'), 'file');
                else
                    mex_exists = exist( fullfile(obj.opts_struct.output_dir,...
                        ['ocp_create.', mexext]), 'file');
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

            % create C object
            try
                obj.C_ocp = ocp_create(obj.model_struct, obj.opts_struct);
            catch ex
                str = sprintf('Exception:\n\t%s\n\t%s\n',ex.identifier,ex.message);
                error(str);
            end

            % generate and compile casadi functions
            if (strcmp(obj.opts_struct.codgen_model, 'true') || strcmp(obj.opts_struct.compile_model, 'true'))
                ocp_generate_casadi_ext_fun(obj.model_struct, obj.opts_struct);
            end

            obj.C_ocp_ext_fun = ocp_create_ext_fun();

            % compile mex with model dependency & set pointers for external functions in model
            obj.C_ocp_ext_fun = ocp_set_ext_fun(obj.C_ocp, obj.C_ocp_ext_fun,...
                                             obj.model_struct, obj.opts_struct);

            % precompute
            ocp_precompute(obj.C_ocp);

            % set parameters to nominal value
            if obj.model_struct.dim_np > 0 && ~isempty(obj.opts_struct.parameter_values)
                obj.set('p', obj.opts_struct.parameter_values)
            end

        end


        function solve(obj)
            ocp_solve(obj.C_ocp);
        end


        function generate_c_code(obj, simulink_opts)
            if nargin < 2
                simulink_opts = get_acados_simulink_opts;
            end
            % set up acados_ocp_nlp_json
            obj.acados_ocp_nlp_json = set_up_acados_ocp_nlp_json(obj, simulink_opts);
            % render templated code
            ocp_generate_c_code(obj);
        end


        function eval_param_sens(obj, field, stage, index)
            ocp_eval_param_sens(obj.C_ocp, field, stage, index);
        end

        function value = get_cost(obj)
            value = ocp_get_cost(obj.C_ocp);
        end

        function set(varargin)
            obj = varargin{1};
            field = varargin{2};
            value = varargin{3};
            if ~isa(field, 'char')
                error('field must be a char vector, use '' ''');
            end
            if nargin==3
                ocp_set(obj.cost_ext_fun_type, obj.cost_ext_fun_type_e, obj.C_ocp, obj.C_ocp_ext_fun, field, value);
            elseif nargin==4
                stage = varargin{4};
                ocp_set(obj.cost_ext_fun_type, obj.cost_ext_fun_type_e, obj.C_ocp, obj.C_ocp_ext_fun, field, value, stage);
            else
                disp('acados_ocp.set: wrong number of input arguments (2 or 3 allowed)');
            end
        end



        function value = get(varargin)
            obj = varargin{1};
            field = varargin{2};
            if ~isa(field, 'char')
                error('field must be a char vector, use '' ''');
            end

            if nargin==2
                value = ocp_get(obj.C_ocp, field);
            elseif nargin==3
                stage = varargin{3};
                value = ocp_get(obj.C_ocp, field, stage);
            else
                disp('acados_ocp.get: wrong number of input arguments (1 or 2 allowed)');
            end
        end

        function [] = store_iterate(varargin)
            %%%  Stores the current iterate of the ocp solver in a json file.
            %%% param1: filename: if not set, use model_name + timestamp + '.json'
            %%% param2: overwrite: if false and filename exists add timestamp to filename

            obj = varargin{1};
            filename = '';
            overwrite = false;

            if nargin>=2
                filename = varargin{2};
                if ~isa(filename, 'char')
                    error('filename must be a char vector, use '' ''');
                end
            end

            if nargin==3
                overwrite = varargin{3};
            end

            if nargin > 3
                disp('acados_ocp.get: wrong number of input arguments (1 or 2 allowed)');
            end

            if strcmp(filename,'')
                filename = [obj.model_struct.name '_iterate.json'];
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
            for i=0:obj.opts_struct.param_scheme_N
                solution.(['x_' num2str(i)]) = obj.get('x', i);
                solution.(['lam_' num2str(i)]) = obj.get('lam', i);
                solution.(['t_' num2str(i)]) = obj.get('t', i);
                solution.(['sl_' num2str(i)]) = obj.get('sl', i);
                solution.(['su_' num2str(i)]) = obj.get('su', i);
            end
            for i=0:obj.opts_struct.param_scheme_N-1
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


        function print(varargin)
            if nargin < 2
                field = 'stat';
            else
                field = varargin{2};
            end

            obj = varargin{1};
            ocp_solver_string = obj.opts_struct.nlp_solver;

            if strcmp(field, 'stat')
                stat = obj.get('stat');
                if strcmp(ocp_solver_string, 'sqp')
                    fprintf('\niter\tres_stat\tres_eq\t\tres_ineq\tres_comp\tqp_stat\tqp_iter\talpha\t');
                    if size(stat,2)>8
                        fprintf('\tqp_res_stat\tqp_res_eq\tqp_res_ineq\tqp_res_comp');
                    end
                    fprintf('\n');
                    for jj=1:size(stat,1)
                        fprintf('%d\t%e\t%e\t%e\t%e\t%d\t%d\t%e', stat(jj,1), stat(jj,2), stat(jj,3), stat(jj,4), stat(jj,5), stat(jj,6), stat(jj,7), stat(jj, 8));
                        if size(stat,2)>8
                            fprintf('\t%e\t%e\t%e\t%e', stat(jj,9), stat(jj,10), stat(jj,11), stat(jj,12));
                        end
                        fprintf('\n');
                    end
                    fprintf('\n');
                elseif strcmp(ocp_solver_string, 'sqp_rti')
                    fprintf('\niter\tqp_status\tqp_iter');
                    if size(stat,2)>3
                        fprintf('\t\tqp_res_stat\tqp_res_eq\tqp_res_ineq\tqp_res_comp');
                    end
                    fprintf('\n');
                    for jj=1:size(stat,1)
                        fprintf('%d\t%d\t\t%d', stat(jj,1), stat(jj,2), stat(jj,3));
                        if size(stat,2)>3
                            fprintf('\t%e\t%e\t%e\t%e', stat(jj,4), stat(jj,5), stat(jj,6), stat(jj,7));
                        end
                        fprintf('\n');
                    end
                end

            else
                fprintf('unsupported field in function print of acados_ocp.print, got %s', field);
                keyboard
            end

        end


        function delete(obj)
            if ~isempty(obj.C_ocp_ext_fun)
                ocp_destroy_ext_fun(obj.model_struct, obj.C_ocp, obj.C_ocp_ext_fun);
            end
            if ~isempty(obj.C_ocp)
                ocp_destroy(obj.C_ocp);
            end
        end


    end % methods

end % class

