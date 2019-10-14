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

function ocp_json = set_up_acados_ocp_nlp_json(obj)
keyboard
    % create
    ocp_json = acados_template_mex.acados_ocp_nlp_json();
    % TODO(andrea): this is temporary. later on the solver_config
    % object will separate from the OCP object

    % general
    ocp_json.solver_config.qp_solver = upper(obj.opts_struct.qp_solver);
    ocp_json.solver_config.integrator_type = upper(obj.opts_struct.sim_method);
    ocp_json.solver_config.nlp_solver_type = upper(obj.opts_struct.nlp_solver);
    ocp_json.dims.N = upper(obj.opts_struct.param_scheme_N);

    ocp_json.solver_config.tf = obj.model_struct.T;
    ocp_json.model.name = obj.model_struct.name;

    %% dims
    % path
    ocp_json.dims.nx = obj.model_struct.dim_nx;
    ocp_json.dims.nu = obj.model_struct.dim_nu;
    ocp_json.dims.nz = obj.model_struct.dim_nz;
    % missing in MEX (?!)
    % ocp_json.dims.np = obj.model_struct.dim_np;
    ocp_json.dims.ny = obj.model_struct.dim_ny;
    ocp_json.dims.nbx = obj.model_struct.dim_nbx;
    ocp_json.dims.nbu = obj.model_struct.dim_nbu;
    ocp_json.dims.ng = obj.model_struct.dim_ng;
    ocp_json.dims.nh = obj.model_struct.dim_nh;
    ocp_json.dims.ns = obj.model_struct.dim_ns;
    ocp_json.dims.nsbx = obj.model_struct.dim_nsbx;
    ocp_json.dims.nsbu = obj.model_struct.dim_nsbu;

    % missing in template
    % ocp_json.dims.nsg = obj.model_struct.nsg;

    % missing in MEX
    % ocp_json.dims.npd = obj.model_struct.npd;
    % ocp_json.dims.npd_e = obj.model_struct.npd_e;

    ocp_json.dims.nsh = obj.model_struct.nsh;

    % terminal
    ocp_json.dims.ng_e = obj.model_struct.ng_e;
    ocp_json.dims.ny_e = obj.model_struct.ny_e;
    ocp_json.dims.nh_e = obj.model_struct.nh_e;
    ocp_json.dims.ns_e = obj.model_struct.ns_e;
    ocp_json.dims.nsh_e = obj.model_struct.nsh_e;
    % missing in MEX
    % ocp_json.dims.nbx_e = obj.model_struct.nbx_e;
    % ocp_json.dims.ns_e = obj.model_struct.ns_e;

    %% types
    ocp_json.cost.cost_type = upper(obj.model_struct.cost_type);
    ocp_json.cost.cost_type_e = upper(obj.model_struct.cost_type_e);
    ocp_json.constraints.constr_type = upper(obj.model_struct.constr_type);
    ocp_json.constraints.constr_type_e = upper(obj.model_struct.constr_type_e);

    %% constraints
    % path
    if isfield(obj.model_struct, 'constr_x0')
        ocp_json.constraints.x0 = obj.model_struct.constr_x0;
    else
        warning('constr_x0 not defined for ocp json.');
    end

    if ocp_json.dims.nbx > 0
        ocp_json.constraints.idxbx = J_to_idx( obj.model_struct.constr_Jbx );
        ocp_json.constraints.lbx = obj.model_struct.constr_lbx;
        ocp_json.constraints.ubx = obj.model_struct.constr_ubx;
    end

    if ocp_json.dims.nbu > 0
        ocp_json.constraints.idxbu = J_to_idx( obj.model_struct.constr_Jbu );
        ocp_json.constraints.lbu = obj.model_struct.constr_lbu;
        ocp_json.constraints.ubu = obj.model_struct.constr_ubu;
    end

    if ocp_json.dims.ng > 0
        ocp_json.constraints.C = obj.model_struct.constr_C;
        ocp_json.constraints.D = obj.model_struct.constr_D;
        ocp_json.constraints.lg = obj.model_struct.constr_lg;
        ocp_json.constraints.ug = obj.model_struct.constr_ug;
    end

    if ocp_json.dims.nh > 0
        ocp_json.con_h.name = 'expr_h';
        ocp_json.con_h.expr = obj.model_struct.constr_expr_h;
        ocp_json.constraints.lh = value;
        ocp_json.constraints.uh = value;
        % TODO(oj): can we get rid of the following?
        if isfield(obj.model_struct, 'sym_x')
            ocp_json.con_h.x = obj.model_struct.sym_x;
        end
        if isfield(obj.model_struct, 'sym_u')
            ocp_json.con_h.u = obj.model_struct.sym_u;
        end
        if isfield(obj.model_struct, 'sym_z')
            ocp_json.con_h.z = obj.model_struct.sym_z;
        end
    end

    if ocp_json.nsbx > 0
        ocp_json.constraints.idxsbx = J_to_idx_slack(obj.model_struct.Jsbx);
        ocp_json.constraints.lsbx = obj.model_struct.lsbx;
        ocp_json.constraints.usbx = obj.model_struct.usbx;
    end

    if ocp_json.nsbu > 0
        ocp_json.constraints.idxsbu = J_to_idx_slack(obj.model_struct.Jsbu);
        ocp_json.constraints.lsbu = obj.model_struct.lsbu;
        ocp_json.constraints.usbu = obj.model_struct.usbu;
    end

    if ocp_json.nsbh > 0
        ocp_json.constraints.idxsh = J_to_idx_slack(obj.model_struct.Jsh);
    end

    if ocp_json.nsg > 0
        error('Jsg not implmented in code-gen backend');
    end


    % terminal
    if ocp_json.dims.ng_e > 0
        ocp_json.constraints.C_e = obj.model_struct.constr_C_e;
        ocp_json.constraints.lg_e = obj.model_struct.constr_lg_e;
        ocp_json.constraints.ug_e = obj.model_struct.constr_ug_e;
    end

    if ocp_json.dims.nh_e > 0    
        ocp_json.con_h_e.name = 'expr_h_e';
        ocp_json.con_h_e.expr = obj.model_struct.constr_expr_h_e;
        ocp_json.constraints.lh_e = obj.model_struct.constr_lh_e;
        ocp_json.constraints.uh_e = obj.model_struct.constr_uh_e;
        % TODO(oj): can we get rid of the following?
        if isfield(obj.model_struct, 'sym_x')
            ocp_json.con_h_e.x = obj.model_struct.sym_x;
        end
        if isfield(obj.model_struct, 'sym_u')
            ocp_json.con_h_e.u = obj.model_struct.sym_u;
        end
        if isfield(obj.model_struct, 'sym_z')
            ocp_json.con_h_e.z = obj.model_struct.sym_z;
        end
    end

    if ocp_json.nsg_e > 0
        error('Jsg_e not implmented in code-gen backend')
    end

    if ocp_json.nsh_e > 0
        ocp_json.constraints.idxsh_e = J_to_idx_slack(obj.model_struct.Jsh_e);
    end

    %% Cost
    ocp_json.cost.Vu = obj.model_struct.cost_Vu;
    ocp_json.cost.Vx = obj.model_struct.cost_Vx;
    ocp_json.cost.Vx_e = obj.model_struct.cost_Vx_e;
    ocp_json.cost.W = obj.model_struct.cost_W;
    ocp_json.cost.Vz = obj.model_struct.cost_Vz;
    ocp_json.cost.W_e = obj.model_struct.cost_W_e;
    ocp_json.cost.yref = obj.model_struct.cost_y_ref;
    ocp_json.cost.yref_e = obj.model_struct.cost_y_ref_e;
    ocp_json.cost.Zl = obj.model_struct.cost_Zl;
    ocp_json.cost.Zl_e = obj.model_struct.cost_Zl_e;
    ocp_json.cost.Zu = obj.model_struct.cost_Zu;
    ocp_json.cost.Zu_e = obj.model_struct.cost_Zu_e;
    ocp_json.cost.zl = obj.model_struct.cost_zl;
    ocp_json.cost.zl_e = obj.model_struct.cost_zl_e;
    ocp_json.cost.zu = obj.model_struct.cost_zu;
    ocp_json.cost.zu_e = obj.model_struct.cost_Zu_e;

    ocp_json.model.x = obj.model_struct.sym_x;
    ocp_json.model.f_impl_expr = value;
    ocp_json.model.f_expl_expr = value;

    if isfield(obj.model_struct, 'sym_u')
        ocp_json.model.u = obj.model_struct.sym_u;
    end
    if isfield(obj.model_struct, 'sym_z')
        ocp_json.model.z = obj.model_struct.sym_z;
    end
    if isfield(obj.model_struct, 'sym_xdot')
        ocp_json.model.xdot = obj.model_struct.sym_xdot;
    end

    keyboard
end


%% auxilary functions

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
        idx(i) = this_idx - 1; % store 0-based index
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