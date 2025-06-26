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

classdef {{ name }}_mex_solver < handle

    properties
        C_ocp
        N
        name
        code_gen_dir
    end % properties



    methods

        % constructor
        function obj = {{ name }}_mex_solver(solver_creation_opts)
            if solver_creation_opts.compile_mex_wrapper
                make_mex_{{ name }}();
            end
            obj.C_ocp = acados_mex_create_{{ name }}();
            % to have path to destructor when changing directory
            addpath('.')
            obj.N = {{ solver_options.N_horizon }};
            obj.name = '{{ name }}';
        end

        % destructor
        function delete(obj)
            disp("delete template...");
            if ~isempty(obj.C_ocp)
                acados_mex_free_{{ name }}(obj.C_ocp);
            end
            disp("done.");
        end

        function solve(obj)
            acados_mex_solve_{{ name }}(obj.C_ocp);
        end

        function status = custom_update(obj, data)
            status = acados_mex_custom_update_{{ name }}(obj.C_ocp, data);
        end

        function set(varargin)
            obj = varargin{1};
            field = varargin{2};
            value = varargin{3};
            if ~isa(field, 'char')
                error('field must be a char vector, use '' ''');
            end
            if nargin==3
                acados_mex_set_{{ name }}(obj.C_ocp, field, value);
            elseif nargin==4
                stage = varargin{4};
                acados_mex_set_{{ name }}(obj.C_ocp, field, value, stage);
            else
                disp('acados_ocp.set: wrong number of input arguments (2 or 3 allowed)');
            end
        end


        function set_params_sparse(varargin)
            % usage:
            % ocp.set_params_sparse(idx_values, param_values, Optional[stage])
            % updates the parameters with indices idx_values (0 based) at stage with the new values new_p_values.
            % if stage is not provided, sparse parameter update is performed for all stages.
            obj = varargin{1};
            idx_values = varargin{2};
            param_values = varargin{3};
            value = [idx_values(:), param_values(:)];
            field = 'params_sparse';

            if ~isa(field, 'char')
                error('field must be a char vector, use '' ''');
            end
            if nargin==3
                acados_mex_set_{{ name }}(obj.C_ocp, field, value);
            elseif nargin==4
                stage = varargin{4};
                acados_mex_set_{{ name }}(obj.C_ocp, field, value, stage);
            else
                disp('acados_ocp.set_params_sparse: wrong number of input arguments (3 or 4 allowed)');
            end
        end

        function value = get_cost(obj)
            value = ocp_get(obj.C_ocp, 'cost_value');
        end

        function value = evaluate_constraints_and_get_violation(obj)
            value = ocp_get(obj.C_ocp, 'constraint_violation');
        end

        function eval_param_sens(obj, field, stage, index)
            ocp_eval_param_sens(obj.C_ocp, field, stage, index);
        end

        function value = get(varargin)
            % usage:
            % obj.get(field, value, [stage])
            obj = varargin{1};
            field = varargin{2};
            if ~isa(field, 'char')
                error('field must be a char vector, use '' ''');
            end
            if any(strfind('sens', field))
                error('field sens* (sensitivities of optimal solution) not yet supported for templated MEX.')
            end

            if nargin==2
                value = ocp_get(obj.C_ocp, field);
            elseif nargin==3
                stage = varargin{3};
                value = ocp_get(obj.C_ocp, field, stage);
            elseif nargin==4
                stage = varargin{3};
                iteration = varargin{4};
                value = ocp_get(obj.C_ocp, field, stage, iteration);
            else
                disp('acados_ocp.get: wrong number of input arguments (1, 2 or 3 allowed)');
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
            for i=0:obj.N
                solution.(['x_' num2str(i)]) = obj.get('x', i);
                solution.(['lam_' num2str(i)]) = obj.get('lam', i);
                solution.(['sl_' num2str(i)]) = obj.get('sl', i);
                solution.(['su_' num2str(i)]) = obj.get('su', i);
            end
            for i=0:obj.N-1
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

        function [] = reset(obj)
            acados_mex_set_{{ name }}(obj.C_ocp, 'reset', 1);
        end


        % print
        function print(varargin)
            if nargin < 2
                field = 'stat';
            else
                field = varargin{2};
            end

            obj = varargin{1};

            if strcmp(field, 'stat')
                stat = obj.get('stat');
                {%- if solver_options.nlp_solver_type == "SQP" %}
                fprintf('\niter\tres_stat\tres_eq\t\tres_ineq\tres_comp\tqp_stat\tqp_iter\talpha');
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
                {%- else %}
                fprintf('\niter\tqp_status\tqp_iter');
                if size(stat,2)>3
                    fprintf('\tqp_res_stat\tqp_res_eq\tqp_res_ineq\tqp_res_comp');
                end
                fprintf('\n');
                for jj=1:size(stat,1)
                    fprintf('%d\t%d\t\t%d', stat(jj,1), stat(jj,2), stat(jj,3));
                    if size(stat,2)>3
                        fprintf('\t%e\t%e\t%e\t%e', stat(jj,4), stat(jj,5), stat(jj,6), stat(jj,7));
                    end
                    fprintf('\n');
                end
                {% endif %}

            else
                fprintf('unsupported field in function print of acados_ocp.print, got %s', field);
                keyboard
            end

        end

    end % methods

end % class

