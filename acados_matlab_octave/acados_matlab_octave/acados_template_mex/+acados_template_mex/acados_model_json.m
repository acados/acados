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

classdef acados_model_json < handle
    properties
        dyn_ext_fun_type
        dyn_generic_source
        dyn_disc_fun_jac_hess
        dyn_disc_fun_jac
        dyn_disc_fun
        dyn_impl_dae_fun
        dyn_impl_dae_fun_jac
        dyn_impl_dae_jac
        f_impl_expr
        f_expl_expr
        f_phi_expr
        x
        xdot
        u
        z
        name
        p
        gnsf
    end
    methods
        function obj = acados_model_json()
            obj.dyn_ext_fun_type = 'casadi';
            obj.dyn_generic_source = [];
            obj.dyn_disc_fun_jac_hess = [];
            obj.dyn_disc_fun_jac = [];
            obj.dyn_disc_fun = [];
            obj.dyn_impl_dae_fun = [];
            obj.dyn_impl_dae_fun_jac = [];
            obj.dyn_impl_dae_jac = [];
            obj.f_impl_expr = [];
            obj.f_expl_expr = [];
            obj.f_phi_expr = [];
            obj.x = [];
            obj.xdot = [];
            obj.u = [];
            obj.z = [];
            obj.name = [];
            obj.p = [];
            obj.gnsf = struct('nontrivial_f_LO', 1, 'purely_linear', 0);
        end
    end
end
