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


classdef GenerateContext < handle
    properties
        p_global
        problem_name
        opts
        casadi_codegen_opts
        list_funname_dir_pairs  % list of (function_name, output_dir) pairs, files that are generated
        generic_funname_dir_pairs % list of (function_name, output_dir) pairs, files that are not generated
        function_input_output_pairs
        function_dyn_cost_constr_types
        casadi_fun_opts
        global_data_sym
        global_data_expr
    end

    methods
        function obj = GenerateContext(p_global, problem_name, opts)
            import casadi.*
            if nargin < 3
                opts = [];
            end
            obj.p_global = p_global;

            if ~(isempty(obj.p_global) || length(obj.p_global) == 0)
                check_casadi_version_supports_p_global();
            end
            obj.problem_name = problem_name;

            obj.global_data_sym = [];
            obj.global_data_expr = [];

            obj.opts = opts;
            obj.casadi_codegen_opts = struct();
            obj.casadi_codegen_opts.mex = false;
            obj.casadi_codegen_opts.casadi_int = 'int';
            obj.casadi_codegen_opts.casadi_real = 'double';
            try
                CodeGenerator('foo', struct('force_canonical', true));
                obj.casadi_codegen_opts.force_canonical = false;
            catch
                % Option does not exist
            end

            obj.list_funname_dir_pairs = {};
            obj.function_input_output_pairs = {};
            obj.function_dyn_cost_constr_types = {};
            obj.generic_funname_dir_pairs = {};

            obj.casadi_fun_opts = struct();

            try
                dummy = MX.sym('dummy');
                cse(dummy); % Check if cse exists
                obj.casadi_fun_opts.cse = true;
            catch
                disp('NOTE: Please consider updating to CasADi 3.6.6 which supports common subexpression elimination. This might speed up external function evaluation.');
            end
        end

        function obj = add_function_definition(obj, name, inputs, outputs, output_dir, dyn_cost_constr_type)
            obj.list_funname_dir_pairs{end+1} = {name, output_dir};
            obj.function_input_output_pairs{end+1} = {inputs, outputs};
            obj.function_dyn_cost_constr_types{end+1} = dyn_cost_constr_type;
        end

        function obj = finalize(obj)
            import casadi.*
            if ~(isempty(obj.p_global) || length(obj.p_global) == 0)
                obj.setup_p_global_precompute_fun()
            end

            obj.generate_functions();
        end


        function obj = add_external_function_file(obj, fun_name, output_dir)
            % remove trailing .c if present
            if endsWith_custom(fun_name, '.c')
                fun_name = fun_name(1:end-2);
            end
            obj.generic_funname_dir_pairs{end+1} = {fun_name, output_dir};
        end

        function out = get_external_function_file_list(obj, ocp_specific)
            out = {};
            for i = 1:numel(obj.generic_funname_dir_pairs)
                fun_name = obj.generic_funname_dir_pairs{i}{1};
                fun_dir = obj.generic_funname_dir_pairs{i}{2};
                rel_fun_dir = relative_path(fun_dir, obj.opts.code_export_directory);
                is_ocp_specific = ~endsWith_custom(rel_fun_dir, 'model');
                if ocp_specific ~= is_ocp_specific
                    continue;
                end
                out{end+1} = [rel_fun_dir, '/', fun_name, '.c'];
            end
            for i = 1:numel(obj.list_funname_dir_pairs)
                fun_name = obj.list_funname_dir_pairs{i}{1};
                fun_dir = obj.list_funname_dir_pairs{i}{2};
                rel_fun_dir = relative_path(fun_dir, obj.opts.code_export_directory);
                is_ocp_specific = ~endsWith_custom(rel_fun_dir, 'model');
                if ocp_specific ~= is_ocp_specific
                    continue;
                end
                out{end+1} = [rel_fun_dir, '/', fun_name, '.c'];
            end
        end

        function n_global_data = get_n_global_data(self)
            n_global_data = length(self.global_data_sym);
        end
    end

    methods (Access = private)

        function [] = setup_p_global_precompute_fun(self)

            import casadi.*

            precompute_pairs = {};

            for i = 1:length(self.function_input_output_pairs)
                outputs = cse(self.function_input_output_pairs{i}{2});

                % detect parametric expressions in p_global
                [outputs_ret, symbols, param_expr] = extract_parametric(outputs, self.p_global, struct('extract_trivial', true));

                % substitute previously detected param_expr in outputs
                symbols_to_add = {};
                param_expr_to_add = {};
                for jj = 1:length(symbols)
                    add = true;
                    sym_new = symbols{jj};
                    expr_new = param_expr{jj};
                    for k = 1:length(precompute_pairs)
                        sym = precompute_pairs{k}{1};
                        expr = precompute_pairs{k}{2};
                        if is_equal(expr, expr_new)
                            for kkk = 1:length(outputs_ret)
                                outputs_ret{kkk} = substitute(outputs_ret{kkk}, sym_new, sym);
                            end
                            add = false;
                            break
                        end
                    end
                    if add
                        symbols_to_add{end+1} = sym_new;
                        param_expr_to_add{end+1} = expr_new;
                    end
                end


                % Update output expressions with the ones that use the extracted expressions
                self.function_input_output_pairs{i}{2} = outputs_ret;

                % Store the new input symbols and extracted expressions
                for j = 1:length(symbols_to_add)
                    precompute_pairs{end+1} = {symbols_to_add{j}, param_expr_to_add{j}};
                end
            end

            % Concatenate global data symbols and expressions
            global_data_sym_list = cellfun(@(pair) pair{1}, precompute_pairs, 'UniformOutput', false);
            self.global_data_sym = vertcat(global_data_sym_list{:});

            global_data_expr_list = cellfun(@(pair) pair{2}, precompute_pairs, 'UniformOutput', false);
            self.global_data_expr = cse(vertcat(global_data_expr_list{:}));

            % make sure global_data is dense
            if length(self.global_data_expr) > 0
                self.global_data_expr = sparsity_cast(self.global_data_expr, Sparsity.dense(self.global_data_expr.nnz()));
                self.global_data_sym = sparsity_cast(self.global_data_sym, Sparsity.dense(self.global_data_sym.nnz()));
            end

            % Assert length match
            assert(length(self.global_data_expr) == length(self.global_data_sym), ...
                   sprintf('Length mismatch: %d != %d', length(self.global_data_expr), length(self.global_data_sym)));

            if length(self.global_data_expr) > 0
                % Add global data as input to all functions
                for i = 1:length(self.function_input_output_pairs)
                    self.function_input_output_pairs{i}{1}{end+1} = self.global_data_sym;
                end

                % Define output directory and function name
                output_dir = fullfile(pwd, self.opts.code_export_directory);
                fun_name = sprintf('%s_p_global_precompute_fun', self.problem_name);

                % Add function definition
                self.add_function_definition(fun_name, {self.p_global}, {self.global_data_expr}, output_dir, 'precompute');
            else
                disp("WARNING: No CasADi function depends on p_global.")
            end


        end


        function [] = generate_functions(obj)
            import casadi.*

            for i = 1:numel(obj.list_funname_dir_pairs)
                name = obj.list_funname_dir_pairs{i}{1};
                output_dir = obj.list_funname_dir_pairs{i}{2};
                inputs = obj.function_input_output_pairs{i}{1};
                outputs = obj.function_input_output_pairs{i}{2};
                dyn_cost_constr_type = obj.function_dyn_cost_constr_types{i};

                % fprintf('Generating function %s in directory %s\n', name, output_dir);
                % disp('Inputs:');
                % disp(inputs);
                % disp('Outputs:');
                % disp(outputs);
                try
                    fun = Function(name, inputs, outputs, obj.casadi_fun_opts);
                catch e
                    fprintf('Error while setting up casadi function %s\n', name);
                    rethrow(e);
                end

                if ((strcmp(dyn_cost_constr_type, 'dyn') && obj.opts.ext_fun_expand_dyn) ...
                    || (strcmp(dyn_cost_constr_type, 'cost') && obj.opts.ext_fun_expand_cost) ...
                    || (strcmp(dyn_cost_constr_type, 'constr') && obj.opts.ext_fun_expand_constr) ...
                    || (strcmp(dyn_cost_constr_type, 'precompute') && obj.opts.ext_fun_expand_precompute))
                    try
                        fun = fun.expand();
                    catch
                        warning(['Failed to expand the CasADi function ' name '.'])
                    end
                end

                % setup and change directory
                cwd = pwd;
                check_dir_and_create(output_dir);
                cd(output_dir);

                % generate function
                try
                    fun.generate(name, obj.casadi_codegen_opts);
                catch e
                    fprintf('Error while generating function %s in directory %s\n', name, output_dir);
                    rethrow(e);
                end

                % change back to original directory
                cd(cwd);
            end
        end
    end
end


function check_casadi_version_supports_p_global()
    import casadi.*
    try
        dummy = MX.sym('dummy');
        % Check if the required functions exist in CasADi
        extract_parametric(dummy, dummy, struct('extract_trivial', true));  % Check if extract_parametric exists
        cse(dummy); % Check if cse exists
        blazing_spline('blazing_spline', {[1, 2, 3], [1, 2, 3]});
    catch
        error('CasADi version does not support extract_parametric or cse functions, thus it is not compatible with p_global in acados. Please install nightly-main release or later, see: https://github.com/casadi/casadi/releases/tag/nightly-main');
    end
end