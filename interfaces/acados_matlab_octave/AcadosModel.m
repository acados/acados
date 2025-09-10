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
        p_global
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

        % GNSF
        sym_gnsf_y
        sym_gnsf_uhat
        dyn_gnsf_A
        dyn_gnsf_A_LO
        dyn_gnsf_B
        dyn_gnsf_B_LO
        dyn_gnsf_E
        dyn_gnsf_E_LO
        dyn_gnsf_C
        dyn_gnsf_c
        dyn_gnsf_c_LO
        dyn_gnsf_L_x
        dyn_gnsf_L_xdot
        dyn_gnsf_L_z
        dyn_gnsf_L_u
        dyn_gnsf_idx_perm_x
        dyn_gnsf_ipiv_x
        dyn_gnsf_idx_perm_z
        dyn_gnsf_ipiv_z
        dyn_gnsf_idx_perm_f
        dyn_gnsf_ipiv_f

        gnsf_nontrivial_f_LO
        gnsf_purely_linear

        dyn_gnsf_expr_phi
        dyn_gnsf_expr_f_lo
    end
    methods
        function obj = AcadosModel()
            obj.name = 'acados_model';
            obj.x = [];
            obj.xdot = [];
            obj.u = [];
            obj.z = [];
            obj.p = [];
            obj.p_global = [];
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

            obj.gnsf_nontrivial_f_LO = 1;
            obj.gnsf_purely_linear = 0;
        end


        function make_consistent(obj, dims)
            import casadi.*
            if isa(obj.x, 'casadi.SX')
                empty_var = SX.sym('empty_var', 0, 0);
                isSX = true;
            elseif isa(obj.x, 'casadi.MX')
                empty_var = MX.sym('empty_var', 0, 0);
                isSX = false;
            else
                error('Unsupported type: model.x must be casadi.SX or casadi.MX');
            end

            if iscolumn(obj.x)
                dims.nx = size(obj.x, 1);
            else
                error('model.x should be column vector of dimension > 0.');
            end

            if isempty(obj.p)
                dims.np = 0;
                obj.p = empty_var;
            elseif iscolumn(obj.p) || (isa(obj.p, 'casadi.SX') == isSX && length(obj.p) == 0)
                dims.np = size(obj.p, 1);
            else
                error('model.p should be column vector.');
            end

            if isempty(obj.p_global)
                dims.np_global = 0;
                obj.p_global = empty_var;
            elseif iscolumn(obj.p_global) || (isa(obj.p_global, 'casadi.SX') == isSX && length(obj.p_global) == 0)
                if any(which_depends(obj.p_global, obj.p))
                    error('model.p_global must not depend on model.p')
                end
                dims.np_global = size(obj.p_global, 1);
            else
                error('model.p_global should be column vector.');
            end

            if isempty(obj.xdot)
                obj.xdot = empty_var;
            elseif ~(isa(obj.xdot, 'casadi.SX') == isSX && length(obj.xdot) == 0) && (~iscolumn(obj.xdot) || size(obj.xdot, 1) ~= dims.nx)
                error('model.xdot should be a column vector of size nx.');
            end

            if isempty(obj.z)
                dims.nz = 0;
                obj.z = empty_var;
            elseif iscolumn(obj.z) || (isa(obj.z, 'casadi.SX') == isSX && length(obj.z) == 0)
                dims.nz = size(obj.z, 1);
            else
                error('model.z should be column vector.');
            end
            if isempty(obj.u)
                dims.nu = 0;
                obj.u = empty_var;
            elseif iscolumn(obj.u) || (isa(obj.u, 'casadi.SX') == isSX && length(obj.u) == 0)
                dims.nu = size(obj.u, 1);
            else
                error('model.u should be column vector.');
            end

            % sanity checks
            vars_and_names = {obj.x, 'x'; obj.xdot, 'xdot'; obj.u, 'u'; obj.z, 'z'; obj.p, 'p'; obj.p_global, 'p_global'};
            for i = 1:size(vars_and_names, 1)
                symbol = vars_and_names{i, 1};
                var_name = vars_and_names{i, 2};
                if ~(isa(symbol, 'casadi.MX') || isa(symbol, 'casadi.SX'))
                    error(['model.' var_name ' must be casadi.MX or casadi.SX, got ' class(symbol)]);
                end
                if ~symbol.is_valid_input()
                    error(['model.' var_name ' must be valid CasADi symbol to be used as input for functions, got ' str(symbol)]);
                end
            end

            % model output dimension nx_next: dimension of the next state
            if isa(dims, 'AcadosOcpDims')
                if ~isempty(obj.disc_dyn_expr)
                    dims.nx_next = length(obj.disc_dyn_expr);
                else
                    dims.nx_next = length(obj.x);
                end
            end

            if ~isempty(obj.f_impl_expr)
                if length(obj.f_impl_expr) ~= (dims.nx + dims.nz)
                    error(sprintf('model.f_impl_expr must have length nx + nz = %d + %d, got %d', dims.nx, dims.nz, length(obj.f_impl_expr)));
                end
            end

            if ~isempty(obj.f_expl_expr)
                if length(obj.f_expl_expr) ~= dims.nx
                    error(sprintf('model.f_expl_expr must have length nx = %d, got %d', dims.nx, length(obj.f_expl_expr)));
                end
            end
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

        function out = convert_to_struct_for_json_dump(self)
            out = struct();
            % all but casadi expressions / variables
            out.name = self.name;
            out.dyn_ext_fun_type = self.dyn_ext_fun_type;
            out.dyn_generic_source = self.dyn_generic_source;
            out.dyn_disc_fun_jac_hess = self.dyn_disc_fun_jac_hess;
            out.dyn_disc_fun_jac = self.dyn_disc_fun_jac;
            out.dyn_disc_fun = self.dyn_disc_fun;
            out.dyn_impl_dae_fun_jac = self.dyn_impl_dae_fun_jac;
            out.dyn_impl_dae_jac = self.dyn_impl_dae_jac;
            out.dyn_impl_dae_fun = self.dyn_impl_dae_fun;

            out.gnsf_nontrivial_f_LO = self.gnsf_nontrivial_f_LO;
            out.gnsf_purely_linear = self.gnsf_purely_linear;
        end
    end
end
