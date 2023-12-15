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

classdef acados_ocp_nlp_json < handle
    properties
        dims
        cost
        constraints
        solver_options
        model
        parameter_values % initial value of the parameter
        acados_include_path
        acados_lib_path
        problem_class
        simulink_opts
        cython_include_dirs
        code_export_directory
        json_file
        shared_lib_ext
        name
    end
    methods
        function obj = acados_ocp_nlp_json(simulink_opts)
            obj.dims = acados_template_mex.ocp_nlp_dims_json();
            obj.cost = acados_template_mex.ocp_nlp_cost_json();
            obj.constraints = acados_template_mex.ocp_nlp_constraints_json();
            obj.solver_options = acados_template_mex.ocp_nlp_solver_options_json();
            obj.model = acados_template_mex.acados_model_json();
            obj.acados_include_path = [];
            obj.acados_lib_path = [];
            obj.parameter_values = [];
            obj.problem_class = 'OCP';
            obj.simulink_opts = simulink_opts;
            obj.cython_include_dirs = [];
            obj.json_file = 'acados_ocp_nlp.json';
            obj.shared_lib_ext = '.so';
            obj.name = 'ocp';
            if ismac()
                obj.shared_lib_ext = '.dylib';
            end
            % obj.code_export_directory;
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
end

