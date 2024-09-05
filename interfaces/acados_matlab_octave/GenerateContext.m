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
        pool_names = {};
        p_global_expressions = {};
        opts
        casadi_codegen_opts = struct();
        list_funname_dir_pairs = cell();  % list of (function_name, output_dir), NOTE: this can be used to simplify template based code generation!
        functions_to_generate = cell();
    end

    methods
        function obj = GenerateContext(p_global, problem_name, opts)
            if nargin < 3
                opts = [];
            end
            obj.p_global = p_global;
            obj.problem_name = problem_name;
            obj.opts = opts;

            obj.casadi_codegen_opts.mex = false;
            obj.casadi_codegen_opts.casadi_int = 'int';
            obj.casadi_codegen_opts.casadi_real = 'double';
        end

        function obj = add_function_definition(obj, name, inputs, outputs, output_dir)
            import casadi.*

            if isempty(obj.p_global) || length(obj.p_global) == 0
                % normal behavior (p_global is empty)
                fun = Function(name, inputs, outputs);
                obj = obj.__add_function(name, output_dir, fun);
            else
                % interesting behavior
                check_casadi_version_supports_p_global();
                inputs_augmented = [inputs, obj.p_global];
                fun = Function(name, inputs_augmented, outputs);

                outputs_call = fun.call(inputs_augmented, struct('always_inline', true, 'never_inline', false));

                [outputs_ret, symbols, param] = extract_parametric(outputs_call, obj.p_global);
                symbols = symbols.primitives();

                pools = {};
                for i = 1:numel(symbols)
                    name_e = [name, '|', num2str(i)];
                    pools{end+1} = MX(DM.zeros(symbols{i}.sparsity()), name_e);
                    obj.pool_names{end+1} = name_e;
                end

                outputs_ret = substitute(outputs_ret, symbols, pools);
                obj.p_global_expressions = [obj.p_global_expressions, param.primitives()];

                fun_mod = Function(fun.name(), inputs, outputs_ret);
                obj = obj.__add_function(name, output_dir, fun_mod);
            end
        end

        function obj = finalize(obj)
            import casadi.*
            if ~(isempty(obj.p_global) || length(obj.p_global) == 0)
                y = cse(obj.p_global_expressions);
                output_dir = obj.opts.code_export_directory;
                fun_name = [obj.problem_name, '_p_global_precompute_fun'];
                fun = Function(fun_name, {obj.p_global}, y, {'p_global'}, obj.pool_names);
                obj = obj.__add_function(fun_name, output_dir, fun);
            end

            obj = obj.__generate_functions();
        end
    end

    methods (Access = private)
        function obj = __add_function(obj, name, output_dir, fun)
            obj.list_funname_dir_pairs{end+1} = {name, output_dir};
            obj.functions_to_generate{end+1} = fun;
        end

        function obj = __generate_functions(obj)

            for i = 1:numel(obj.list_funname_dir_pairs)
                name = obj.list_funname_dir_pairs{i}{1};
                output_dir = obj.list_funname_dir_pairs{i}{2};
                fun = obj.functions_to_generate{i};

                % setup and change directory
                cwd = pwd;
                if ~exist(output_dir, 'dir')
                    mkdir(output_dir);
                end
                cd(output_dir);

                % Generate function
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
        % Check if the required functions exist in CasADi
        extract_parametric;  % Check if extract_parametric exists
        cse;                 % Check if cse exists
    catch
        error('CasADi version does not support extract_parametric or cse functions.\nNeeds nightly-se release or later, see: https://github.com/casadi/casadi/releases/tag/nightly-se');
    end
end