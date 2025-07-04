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



classdef AcadosSimOptions < handle
    properties
        integrator_type
        collocation_type
        Tsim
        num_stages
        num_steps
        newton_iter
        newton_tol
        jac_reuse
        sens_forw
        sens_adj
        sens_algebraic
        sens_hess
        output_z
        ext_fun_compile_flags
        ext_fun_expand_dyn
        compile_interface
        with_batch_functionality
    end

    methods
        function obj = AcadosSimOptions()
            obj.integrator_type = 'ERK';
            obj.collocation_type = 'GAUSS_LEGENDRE';
            obj.Tsim = [];
            obj.num_stages = 4;
            obj.num_steps = 1;
            obj.newton_iter = 3;
            obj.newton_tol = 0.;
            obj.sens_forw = true;
            obj.sens_adj = false;
            obj.sens_algebraic = false;
            obj.sens_hess = false;
            obj.output_z = true;
            obj.jac_reuse = 0;
            % check whether flags are provided by environment variable
            env_var = getenv("ACADOS_EXT_FUN_COMPILE_FLAGS");
            if isempty(env_var)
                obj.ext_fun_compile_flags = '-O2';
            else
                obj.ext_fun_compile_flags = env_var;
            end
            obj.ext_fun_expand_dyn = false;
            obj.with_batch_functionality = false;
            obj.compile_interface = []; % corresponds to automatic detection, possible values: true, false, []
        end

        function s = struct(self)
            if exist('properties')
                publicProperties = eval('properties(self)');
            else
                publicProperties = fieldnames(self);
            end
            s = struct();
            for fi = 1:numel(publicProperties)
                property_name = publicProperties{fi};
                if strcmp(property_name, 'num_stages') || strcmp(property_name, 'num_steps') || strcmp(property_name, 'newton_iter') || ...
                     strcmp(property_name, 'jac_reuse') || strcmp(property_name, 'newton_tol')
                    out_name = strcat('sim_method_', property_name);
                    s.(out_name) = self.(property_name);
                else
                    s.(property_name) = self.(property_name);
                end
            end
        end
    end
end
