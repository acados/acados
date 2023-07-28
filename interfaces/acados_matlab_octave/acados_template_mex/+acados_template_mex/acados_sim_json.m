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

classdef acados_sim_json < handle

    properties
        % file structure
        acados_include_path
        acados_lib_path
        shared_lib_ext
        json_file
        cython_include_dirs
        code_export_directory
        % struct / object
        dims
        model
        sim_options
        % plain data
        parameter_values
        problem_class
    end

    methods
        function obj = acados_sim_json()
            % most fields are initialized as a placeholder
            obj.acados_include_path = [];
            obj.acados_lib_path = [];
            obj.shared_lib_ext = [];
            obj.json_file = [];
            obj.cython_include_dirs = [];
            obj.code_export_directory = [];

            obj.dims = struct( ...
                'nx', [], ...
                'nu', [], ...
                'nz', [], ...
                'np', [] ...
            );
            obj.model = acados_template_mex.acados_model_json();
            obj.sim_options = struct( ...
                'integrator_type', [], ...
                'collocation_type', [], ...
                'sim_method_num_stages', [], ...
                'sim_method_num_steps', [], ...
                'sim_method_newton_iter', [], ...
                'sim_method_newton_tol', [], ...
                'Tsim', [], ...
                'sens_forw', [], ...
                'sens_adj', [], ...
                'sens_algebraic', [], ...
                'sens_hess', [], ...
                'output_z', [], ...
                'sim_method_jac_reuse', [], ...
                'ext_fun_compile_flags', [] ...
            );

            obj.parameter_values = [];
            obj.problem_class = 'SIM';
        end

        function s = struct(self)
            if exist('properties')
                publicProperties = eval('properties(self)');
            else
                publicProperties = fieldnames(self);
            end
            s = struct();
            for fi = 1:numel(publicProperties)
                s.(publicProperties{fi}) = self.(publicProperties{fi});
            end
        end
    end

end % class
