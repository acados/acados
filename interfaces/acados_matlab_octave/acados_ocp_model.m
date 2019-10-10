%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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

classdef acados_ocp_model < handle

    properties
        model_struct
        acados_ocp_nlp_json
    end % properties


    methods(Static)
        function idx = J_to_idx(J)
            size_J = size(J);
            nrows = size_J(1);
            idx = zeros(nrows,1);
            for i = 1:nrows
                this_idx = find(J(i,:));
                if length(this_idx) ~= 1
                    error(['J_to_idx: Invalid J matrix. Exiting. Found more than one nonzero in row ' num2str(i)]);
                end
                if J(i,this_idx) ~= 1
                    error(['J_to_idx: J matrices can only contain 1s, got J(' num2str(i) ', ' num2str(this_idx) ') = ' num2str(J(i,this_idx)) ]);
                end
                idx(i) = this_idx - 1; % strore 0-based index
            end
        end

        function idx = J_to_idx_slack(J)
            size_J = size(J);
            nrows = size_J(1);
            idx = zeros(nrows,1);
            i_idx = 1;
            for i = 1:nrows
                this_idx = find(J(i,:));
                if length(this_idx) == 1
                    i_idx = i_idx + 1;
                    idx(i_idx) = this_idx - 1; % strore 0-based index
                elseif length(this_idx) > 1
                    error(['J_to_idx: Invalid J matrix. Exiting. Found more than one nonzero in row ' num2str(i)]);
                end
                if J(i,this_idx) ~= 1
                    error(['J_to_idx: J matrices can only contain 1s, got J(' num2str(i) ', ' num2str(this_idx) ') = ' num2str(J(i,this_idx)) ]);
                end
            end
        end
    end % methods(static)
    
    methods
        

        function obj = acados_ocp_model()
            % model structure
            obj.model_struct = struct;
            % default values
            obj.model_struct.name = 'ocp_model';
            obj.model_struct.cost_type = 'linear_ls';
            obj.model_struct.cost_type_e = 'linear_ls';
            obj.model_struct.dyn_type = 'implicit';
            obj.model_struct.constr_type = 'bgh';
            % JSON data
            obj.acados_ocp_nlp_json = acados_template_mex.acados_ocp_nlp_json();
        end



        function obj = set(obj, field, value)

            % check for module name
            tokens = strsplit(field, '_');

            % symbolics
            if (strcmp(tokens{1}, 'sym'))

                if (strcmp(field, 'sym_x'))
                    obj.model_struct.sym_x = value;
                elseif (strcmp(field, 'sym_xdot'))
                    obj.model_struct.sym_xdot = value;
                elseif (strcmp(field, 'sym_u'))
                    obj.model_struct.sym_u = value;
                elseif (strcmp(field, 'sym_z'))
                    obj.model_struct.sym_z = value;
                elseif (strcmp(field, 'sym_p'))
                    obj.model_struct.sym_p = value;
                else
                    disp(['acados_ocp_model: set: wrong field: ', field]);
                    keyboard;
                end

            % cost
            elseif (strcmp(tokens{1}, 'cost'))

                if (strcmp(field, 'cost_type'))
                    obj.model_struct.cost_type = value;
                    obj.acados_ocp_nlp_json.cost.cost_type = upper(value);
                elseif (strcmp(field, 'cost_type_e'))
                    obj.model_struct.cost_type_e = value;
                    obj.acados_ocp_nlp_json.cost.cost_type_e = upper(value);
                elseif (strcmp(field, 'cost_expr_y'))
                    obj.model_struct.cost_expr_y = value;
                elseif (strcmp(field, 'cost_expr_y_e'))
                    obj.model_struct.cost_expr_y_e = value;
                elseif (strcmp(field, 'cost_expr_ext_cost'))
                    obj.model_struct.cost_expr_ext_cost = value;
                elseif (strcmp(field, 'cost_expr_ext_cost_e'))
                    obj.model_struct.cost_expr_ext_cost_e = value;
                elseif (strcmp(field, 'cost_Vu'))
                    obj.model_struct.cost_Vu = value;
                    obj.acados_ocp_nlp_json.cost.Vu = value;
                elseif (strcmp(field, 'cost_Vx'))
                    obj.model_struct.cost_Vx = value;
                    obj.acados_ocp_nlp_json.cost.Vx = value;
                elseif (strcmp(field, 'cost_Vx_e'))
                    obj.model_struct.cost_Vx_e = value;
                    obj.acados_ocp_nlp_json.cost.Vx_e = value;
                elseif (strcmp(field, 'cost_Vz'))
                    obj.model_struct.cost_Vz = value;
                    obj.acados_ocp_nlp_json.cost.Vz = value;
                elseif (strcmp(field, 'cost_W'))
                    obj.model_struct.cost_W = value;
                    obj.acados_ocp_nlp_json.cost.W = value;
                elseif (strcmp(field, 'cost_W_e'))
                    obj.model_struct.cost_W_e = value;
                    obj.acados_ocp_nlp_json.cost.W_e = value;
                elseif (strcmp(field, 'cost_y_ref'))
                    obj.model_struct.cost_y_ref = value;
                    obj.acados_ocp_nlp_json.cost.yref = value;
                elseif (strcmp(field, 'cost_y_ref_e'))
                    obj.model_struct.cost_y_ref_e = value;
                    obj.acados_ocp_nlp_json.cost.yref_e = value;
                elseif (strcmp(field, 'cost_Z'))
                    obj.model_struct.cost_Z = value;
                elseif (strcmp(field, 'cost_Z_e'))
                    obj.model_struct.cost_Z_e = value;
                elseif (strcmp(field, 'cost_Zl'))
                    obj.model_struct.cost_Zl = value;
                    obj.acados_ocp_nlp_json.cost.Zl = value;
                elseif (strcmp(field, 'cost_Zl_e'))
                    obj.model_struct.cost_Zl_e = value;
                    obj.acados_ocp_nlp_json.cost.Zl_e = value;
                elseif (strcmp(field, 'cost_Zu'))
                    obj.model_struct.cost_Zu = value;
                    obj.acados_ocp_nlp_json.cost.Zu = value;
                elseif (strcmp(field, 'cost_Zu_e'))
                    obj.model_struct.cost_Zu_e = value;
                    obj.acados_ocp_nlp_json.cost.Zu_e = value;
                elseif (strcmp(field, 'cost_zl'))
                    obj.model_struct.cost_zl = value;
                    obj.acados_ocp_nlp_json.cost.zl = value;
                elseif (strcmp(field, 'cost_zl_e'))
                    obj.model_struct.cost_zl_e = value;
                    obj.acados_ocp_nlp_json.cost.zl_e = value;
                elseif (strcmp(field, 'cost_z'))
                    obj.model_struct.cost_z = value;
                elseif (strcmp(field, 'cost_z_e'))
                    obj.model_struct.cost_z_e = value;
                elseif (strcmp(field, 'cost_zu'))
                    obj.model_struct.cost_zu = value;
                    obj.acados_ocp_nlp_json.cost.zu = value;
                elseif (strcmp(field, 'cost_zu_e'))
                    obj.model_struct.cost_zu_e = value;
                    obj.acados_ocp_nlp_json.cost.zu_e = value;
                else
                    disp(['acados_ocp_model: set: wrong field: ', field]);
                    keyboard;
                end

            % dynamics
            elseif (strcmp(tokens{1}, 'dyn'))

                if (strcmp(field, 'dyn_type'))
                    obj.model_struct.dyn_type = value;
                elseif (strcmp(field, 'dyn_expr_f'))
                    obj.model_struct.dyn_expr_f = value;
                    if isfield(obj.model_struct, 'sym_x')
                        obj.acados_ocp_nlp_json.model.x = obj.model_struct.sym_x;
                    end
                    if isfield(obj.model_struct, 'sym_u')
                        obj.acados_ocp_nlp_json.model.u = obj.model_struct.sym_u;
                    end
                    if isfield(obj.model_struct, 'sym_z')
                        obj.acados_ocp_nlp_json.model.z = obj.model_struct.sym_z;
                    end
                    if isfield(obj.model_struct, 'sym_xdot')
                        obj.acados_ocp_nlp_json.model.xdot = obj.model_struct.sym_xdot;
                    end
                    % obj.acados_ocp_nlp_json.model.name = 'dyn';
                    obj.acados_ocp_nlp_json.model.f_expl_expr = value;
                elseif (strcmp(field, 'dyn_expr_phi'))
                    obj.model_struct.dyn_expr_phi = value;
                    if isfield(obj.model_struct, 'sym_x')
                        obj.acados_ocp_nlp_json.model.x = obj.model_struct.sym_x;
                    end
                    if isfield(obj.model_struct, 'sym_u')
                        obj.acados_ocp_nlp_json.model.u = obj.model_struct.sym_u;
                    end
                    if isfield(obj.model_struct, 'sym_z')
                        obj.acados_ocp_nlp_json.model.z = obj.model_struct.sym_z;
                    end
                    if isfield(obj.model_struct, 'sym_xdot')
                        obj.acados_ocp_nlp_json.model.xdot = obj.model_struct.sym_xdot;
                    end
                    % obj.acados_ocp_nlp_json.model.name = 'dyn';
                    obj.acados_ocp_nlp_json.model.f_impl_expr = value;
                else
                    disp(['acados_ocp_model: set: wrong field: ', field]);
                    keyboard;
                end

            % constraints
            elseif (strcmp(tokens{1}, 'constr'))

                if (strcmp(field, 'constr_type'))
                    obj.model_struct.constr_type = value;
                    obj.acados_ocp_nlp_json.constraints.constr_type = upper(value);
                elseif (strcmp(field, 'constr_x0'))
                    obj.model_struct.constr_x0 = value;
                    obj.acados_ocp_nlp_json.constraints.x0 = value;
                elseif (strcmp(field, 'constr_Jbx'))
                    obj.model_struct.constr_Jbx = value;
                       obj.acados_ocp_nlp_json.constraints.idxbx = obj.J_to_idx(value);
                elseif (strcmp(field, 'constr_lbx'))
                    obj.model_struct.constr_lbx = value;
                    obj.acados_ocp_nlp_json.constraints.lbx = value;
                elseif (strcmp(field, 'constr_ubx'))
                    obj.model_struct.constr_ubx = value;
                    obj.acados_ocp_nlp_json.constraints.ubx = value;
                elseif (strcmp(field, 'constr_Jbu'))
                    obj.model_struct.constr_Jbu = value;
                      obj.acados_ocp_nlp_json.constraints.idxbu = obj.J_to_idx(value);
                elseif (strcmp(field, 'constr_lbu'))
                    obj.model_struct.constr_lbu = value;
                    obj.acados_ocp_nlp_json.constraints.lbu = value;
                elseif (strcmp(field, 'constr_ubu'))
                    obj.model_struct.constr_ubu = value;
                    obj.acados_ocp_nlp_json.constraints.ubu = value;
                elseif (strcmp(field, 'constr_C'))
                    obj.model_struct.constr_C = value;
                    obj.acados_ocp_nlp_json.constraints.C = value;
                elseif (strcmp(field, 'constr_D'))
                    obj.model_struct.constr_D = value;
                    obj.acados_ocp_nlp_json.constraints.D = value;
                elseif (strcmp(field, 'constr_lg'))
                    obj.model_struct.constr_lg = value;
                    obj.acados_ocp_nlp_json.constraints.lg = value;
                elseif (strcmp(field, 'constr_ug'))
                    obj.model_struct.constr_ug = value;
                    obj.acados_ocp_nlp_json.constraints.ug = value;
                elseif (strcmp(field, 'constr_C_e'))
                    obj.model_struct.constr_C_e = value;
                    obj.acados_ocp_nlp_json.constraints.C_e = value;
                elseif (strcmp(field, 'constr_lg_e'))
                    obj.model_struct.constr_lg_e = value;
                    obj.acados_ocp_nlp_json.constraints.lg_e = value;
                elseif (strcmp(field, 'constr_ug_e'))
                    obj.model_struct.constr_ug_e = value;
                    obj.acados_ocp_nlp_json.constraints.ug_e = value;
                elseif (strcmp(field, 'constr_expr_h'))
                    obj.model_struct.constr_expr_h = value;
                    if isfield(obj.model_struct, 'sym_x')
                        obj.acados_ocp_nlp_json.con_h.x = obj.model_struct.sym_x;
                    end
                    if isfield(obj.model_struct, 'sym_u')
                        obj.acados_ocp_nlp_json.con_h.u = obj.model_struct.sym_u;
                    end
                    if isfield(obj.model_struct, 'sym_z')
                        obj.acados_ocp_nlp_json.con_h.z = obj.model_struct.sym_z;
                    end
                    obj.acados_ocp_nlp_json.con_h.name = 'expr_h';
                    obj.acados_ocp_nlp_json.con_h.expr = value;
                elseif (strcmp(field, 'constr_lh'))
                    obj.model_struct.constr_lh = value;
                    obj.acados_ocp_nlp_json.constraints.lh = value;
                elseif (strcmp(field, 'constr_uh'))
                    obj.model_struct.constr_uh = value;
                    obj.acados_ocp_nlp_json.constraints.uh = value;
                elseif (strcmp(field, 'constr_expr_h_e'))
                    obj.model_struct.constr_expr_h_e = value;
                    if isfield(obj.model_struct, 'sym_x')
                        obj.acados_ocp_nlp_json.con_h_e.x = obj.model_struct.sym_x;
                    end
                    if isfield(obj.model_struct, 'sym_u')
                        obj.acados_ocp_nlp_json.con_h_e.u = obj.model_struct.sym_u;
                    end
                    if isfield(obj.model_struct, 'sym_z')
                        obj.acados_ocp_nlp_json.con_h_e.z = obj.model_struct.sym_z;
                    end
                    obj.acados_ocp_nlp_json.con_h_e.name = 'expr_h_e';
                    obj.acados_ocp_nlp_json.con_h_e.expr = value;
                elseif (strcmp(field, 'constr_lh_e'))
                    obj.model_struct.constr_lh_e = value;
                    obj.acados_ocp_nlp_json.constraints.lh_e = value;
                elseif (strcmp(field, 'constr_uh_e'))
                    obj.model_struct.constr_uh_e = value;
                    obj.acados_ocp_nlp_json.constraints.uh_e = value;
                elseif (strcmp(field, 'constr_Jsbu'))
                    obj.model_struct.constr_Jsbu = value;
                       obj.acados_ocp_nlp_json.constraints.idxsbu = obj.J_to_idx_slack(value);
    %            elseif (strcmp(field, 'constr_lsbu'))
    %                obj.model_struct.constr_lsbu = value;
    %            elseif (strcmp(field, 'constr_usbu'))
    %                obj.model_struct.constr_usbu = value;
                elseif (strcmp(field, 'constr_Jsbx'))
                    obj.model_struct.constr_Jsbx = value;
                    obj.acados_ocp_nlp_json.constraints.idxsbx = obj.J_to_idx_slack(value);
    %            elseif (strcmp(field, 'constr_lsbx'))
    %                obj.model_struct.constr_lsbx = value;
    %            elseif (strcmp(field, 'constr_usbx'))
    %                obj.model_struct.constr_usbx = value;
                elseif (strcmp(field, 'constr_Jsg'))
                    obj.model_struct.constr_Jsg = value;
                    warning('Jsg not implmented in code-gen backend')
    %            elseif (strcmp(field, 'constr_lsg'))
    %                obj.model_struct.constr_lsg = value;
    %            elseif (strcmp(field, 'constr_usg'))
    %                obj.model_struct.constr_usg = value;
                elseif (strcmp(field, 'constr_Jsg_e'))
                    obj.model_struct.constr_Jsg_e = value;
                    warning('Jsg not implmented in code-gen backend')
    %            elseif (strcmp(field, 'constr_lsg_e'))
    %                obj.model_struct.constr_lsg_e = value;
    %            elseif (strcmp(field, 'constr_usg_e'))
    %                obj.model_struct.constr_usg_e = value;
                elseif (strcmp(field, 'constr_Jsh'))
                    obj.model_struct.constr_Jsh = value;
                    obj.acados_ocp_nlp_json.constraints.idxsh = obj.J_to_idx_slack(value);
    %            elseif (strcmp(field, 'constr_lsh'))
    %                obj.model_struct.constr_lsh = value;
    %            elseif (strcmp(field, 'constr_ush'))
    %                obj.model_struct.constr_ush = value;
                elseif (strcmp(field, 'constr_Jsh_e'))
                    obj.model_struct.constr_Jsh_e = value;
                    obj.acados_ocp_nlp_json.constraints.idxsh_e = obj.J_to_idx_slack(value);
    %            elseif (strcmp(field, 'constr_lsh_e'))
    %                obj.model_struct.constr_lsh_e = value;
    %            elseif (strcmp(field, 'constr_ush_e'))
    %                obj.model_struct.constr_ush_e = value;
                else
                    disp(['acados_ocp_model: set: wrong field: ', field]);
                    keyboard;
                end

            % dims
            elseif (strcmp(tokens{1}, 'dim'))

                if (strcmp(field, 'dim_nx'))
                    obj.model_struct.dim_nx = value;
                    obj.acados_ocp_nlp_json.dims.nx = value;
                elseif (strcmp(field, 'dim_nu'))
                    obj.model_struct.dim_nu = value;
                    obj.acados_ocp_nlp_json.dims.nu = value;
                elseif (strcmp(field, 'dim_nz'))
                    obj.model_struct.dim_nz = value;
                    obj.acados_ocp_nlp_json.dims.nz = value;
                elseif (strcmp(field, 'dim_ny'))
                    obj.model_struct.dim_ny = value;
                    obj.acados_ocp_nlp_json.dims.ny = value;
                elseif (strcmp(field, 'dim_ny_e'))
                    obj.model_struct.dim_ny_e = value;
                    obj.acados_ocp_nlp_json.dims.ny_e = value;
                elseif (strcmp(field, 'dim_nbx'))
                    obj.model_struct.dim_nbx = value;
                    obj.acados_ocp_nlp_json.dims.nbx = value;
                elseif (strcmp(field, 'dim_nbu'))
                    obj.model_struct.dim_nbu = value;
                    obj.acados_ocp_nlp_json.dims.nbu = value;
                elseif (strcmp(field, 'dim_ng'))
                    obj.model_struct.dim_ng = value;
                    obj.acados_ocp_nlp_json.dims.ng = value;
                elseif (strcmp(field, 'dim_ng_e'))
                    obj.model_struct.dim_ng_e = value;
                    obj.acados_ocp_nlp_json.dims.ng_e = value;
                elseif (strcmp(field, 'dim_nh'))
                    obj.model_struct.dim_nh = value;
                    obj.acados_ocp_nlp_json.dims.nh = value;
                elseif (strcmp(field, 'dim_nh_e'))
                    obj.model_struct.dim_nh_e = value;
                    obj.acados_ocp_nlp_json.dims.nh_e = value;
                elseif (strcmp(field, 'dim_ns'))
                    obj.model_struct.dim_ns = value;
                    obj.acados_ocp_nlp_json.dims.ns = value;
                elseif (strcmp(field, 'dim_ns_e'))
                    obj.model_struct.dim_ns_e = value;
                    obj.acados_ocp_nlp_json.dims.ns_e = value;
                elseif (strcmp(field, 'dim_nsbu'))
                    obj.model_struct.dim_nsbu = value;
                elseif (strcmp(field, 'dim_nsbx'))
                    obj.model_struct.dim_nsbx = value;
                elseif (strcmp(field, 'dim_nsg'))
                    obj.model_struct.dim_nsg = value;
                elseif (strcmp(field, 'dim_nsg_e'))
                    obj.model_struct.dim_nsg_e = value;
                elseif (strcmp(field, 'dim_nsh'))
                    obj.model_struct.dim_nsh = value;
                    obj.acados_ocp_nlp_json.dims.nsh = value;
                elseif (strcmp(field, 'dim_nsh_e'))
                    obj.model_struct.dim_nsh_e = value;
                    obj.acados_ocp_nlp_json.dims.nsh_e = value;
                elseif (strcmp(field, 'dim_np'))
                    obj.model_struct.dim_np = value;
                else
                    disp(['acados_ocp_model: set: wrong field: ', field]);
                    keyboard;
                end

            % others
            else

                if (strcmp(field, 'name'))
                    obj.model_struct.name = value;
                    obj.acados_ocp_nlp_json.model.name = value;
                elseif (strcmp(field, 'T'))
                    obj.model_struct.T = value;
                    obj.acados_ocp_nlp_json.solver_config.tf = value;
                else
                    disp(['acados_ocp_model: set: wrong field: ', field]);
                    keyboard;
                end
            end    
        end

    end % methods
    




end % class


