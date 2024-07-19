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

classdef AcadosModel < handle
    properties
        name
        x
        xdot
        u
        z
        p
        t
        f_impl_expr
        f_expl_expr
        disc_dyn_expr
        dyn_ext_fun_type
        dyn_generic_source
        dyn_disc_fun_jac_hess
        dyn_disc_fun_jac
        dyn_disc_fun
        dyn_impl_dae_fun_jac
        dyn_impl_dae_jac
        dyn_impl_dae_fun
        gnsf

        con_h_expr_0
        con_phi_expr_0
        con_r_expr_0
        con_r_in_phi_0

        con_h_expr
        con_phi_expr
        con_r_expr
        con_r_in_phi

        con_h_expr_e
        con_phi_expr_e
        con_r_expr_e
        con_r_in_phi_e

        cost_y_expr_0
        cost_y_expr
        cost_y_expr_e

        cost_expr_ext_cost_0
        cost_expr_ext_cost
        cost_expr_ext_cost_e

        cost_expr_ext_cost_custom_hess_0
        cost_expr_ext_cost_custom_hess
        cost_expr_ext_cost_custom_hess_e

        cost_psi_expr_0
        cost_psi_expr
        cost_psi_expr_e

        cost_r_in_psi_expr_0
        cost_r_in_psi_expr
        cost_r_in_psi_expr_e

        cost_conl_custom_outer_hess_0
        cost_conl_custom_outer_hess
        cost_conl_custom_outer_hess_e
    end
    methods
        function obj = AcadosModel()
            obj.name = [];
            obj.x = [];
            obj.xdot = [];
            obj.u = [];
            obj.z = [];
            obj.p = [];
            obj.t = [];

            obj.f_impl_expr = [];
            obj.f_expl_expr = [];
            obj.disc_dyn_expr = [];

            obj.dyn_ext_fun_type = 'casadi';
            obj.dyn_generic_source = [];
            obj.dyn_disc_fun_jac_hess = [];
            obj.dyn_disc_fun_jac = [];
            obj.dyn_disc_fun = [];
            obj.dyn_impl_dae_fun_jac = [];
            obj.dyn_impl_dae_jac = [];
            obj.dyn_impl_dae_fun = [];

            obj.gnsf = struct('nontrivial_f_LO', 1, 'purely_linear', 0);

            obj.con_h_expr_0 = [];
            obj.con_phi_expr_0 = [];
            obj.con_r_expr_0 = [];
            obj.con_r_in_phi_0 = [];

            obj.con_h_expr = [];
            obj.con_phi_expr = [];
            obj.con_r_expr = [];
            obj.con_r_in_phi = [];

            obj.con_h_expr_e = [];
            obj.con_phi_expr_e = [];
            obj.con_r_expr_e = [];
            obj.con_r_in_phi_e = [];

            obj.cost_y_expr_0 = [];
            obj.cost_y_expr = [];
            obj.cost_y_expr_e = [];

            obj.cost_expr_ext_cost_0 = [];
            obj.cost_expr_ext_cost = [];
            obj.cost_expr_ext_cost_e = [];

            obj.cost_expr_ext_cost_custom_hess_0 = [];
            obj.cost_expr_ext_cost_custom_hess = [];
            obj.cost_expr_ext_cost_custom_hess_e = [];

            obj.cost_psi_expr_0 = [];
            obj.cost_psi_expr = [];
            obj.cost_psi_expr_e = [];

            obj.cost_r_in_psi_expr_0 = [];
            obj.cost_r_in_psi_expr = [];
            obj.cost_r_in_psi_expr_e = [];

            obj.cost_conl_custom_outer_hess_0 = [];
            obj.cost_conl_custom_outer_hess = [];
            obj.cost_conl_custom_outer_hess_e = [];
        end
    end
end
