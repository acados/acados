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
        num_threads_in_batch_solve
        compile_interface
    end

    properties (Access = private)
        sim_method_num_stages
        sim_method_num_steps
        sim_method_newton_iter
        sim_method_newton_tol
        sim_method_jac_reuse
    end
    methods
        %% get and set methods
        % num_stages
        function obj = set.num_stages(obj, val)
            obj.sim_method_num_stages = val;
        end
        function val = get.num_stages(obj)
            val = obj.sim_method_num_stages;
        end
        % num_steps
        function obj = set.num_steps(obj, val)
            obj.sim_method_num_steps = val;
        end
        function val = get.num_steps(obj)
            val = obj.sim_method_num_steps;
        end
        % newton_iter
        function obj = set.newton_iter(obj, val)
            obj.sim_method_newton_iter = val;
        end
        function val = get.newton_iter(obj)
            val = obj.sim_method_newton_iter;
        end
        % newton_tol
        function obj = set.newton_tol(obj, val)
            obj.sim_method_newton_tol = val;
        end
        function val = get.newton_tol(obj)
            val = obj.sim_method_newton_tol;
        end
        % jac_reuse
        function obj = set.jac_reuse(obj, val)
            obj.sim_method_jac_reuse = val;
        end
        function val = get.jac_reuse(obj)
            val = obj.sim_method_jac_reuse;
        end
        function obj = AcadosSimOptions()
            obj.integrator_type = 'ERK';
            obj.collocation_type = 'GAUSS_LEGENDRE';
            obj.Tsim = [];
            obj.sim_method_num_stages = 4;
            obj.sim_method_num_steps = 1;
            obj.sim_method_newton_iter = 3;
            obj.sim_method_newton_tol = 0.;
            obj.sens_forw = true;
            obj.sens_adj = false;
            obj.sens_algebraic = false;
            obj.sens_hess = false;
            obj.output_z = true;
            obj.sim_method_jac_reuse = 0;
            obj.ext_fun_compile_flags = '-O2';
            obj.num_threads_in_batch_solve = 1;
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
                s.(publicProperties{fi}) = self.(publicProperties{fi});
            end
        end
    end
end
