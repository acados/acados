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

    model = obj.model_struct;
    % create
    ocp_json = acados_template_mex.acados_ocp_nlp_json();
    % TODO(andrea): this is temporary. later on the solver_options
    % object will separate from the OCP object

    % general
    ocp_json.dims.N = obj.opts_struct.param_scheme_N;
    ocp_json.solver_options.tf = model.T;
    ocp_json.model.name = model.name;
    % modules
    ocp_json.solver_options.qp_solver = upper(obj.opts_struct.qp_solver);
    ocp_json.solver_options.integrator_type = upper(obj.opts_struct.sim_method);
    ocp_json.solver_options.nlp_solver_type = upper(obj.opts_struct.nlp_solver);
    % options
    ocp_json.solver_options.sim_method_num_steps = obj.opts_struct.sim_method_num_steps;
    ocp_json.solver_options.sim_method_num_stages = obj.opts_struct.sim_method_num_stages;
    ocp_json.solver_options.sim_method_newton_iter = obj.opts_struct.sim_method_newton_iter;
    ocp_json.solver_options.nlp_solver_max_iter = obj.opts_struct.nlp_solver_max_iter;
    ocp_json.solver_options.nlp_solver_tol_stat = obj.opts_struct.nlp_solver_tol_stat;
    ocp_json.solver_options.nlp_solver_tol_eq = obj.opts_struct.nlp_solver_tol_eq;
    ocp_json.solver_options.nlp_solver_tol_ineq = obj.opts_struct.nlp_solver_tol_ineq;
    ocp_json.solver_options.nlp_solver_tol_comp = obj.opts_struct.nlp_solver_tol_comp;
    ocp_json.solver_options.nlp_solver_step_length = obj.opts_struct.nlp_solver_step_length;


    %% dims
    nx = model.dim_nx;
    % path
    ocp_json.dims.nx = model.dim_nx;
    ocp_json.dims.nu = model.dim_nu;
    ocp_json.dims.nz = model.dim_nz;
    ocp_json.dims.np = model.dim_np;
    ocp_json.dims.ny = model.dim_ny;
    ocp_json.dims.nbx = model.dim_nbx;
    ocp_json.dims.nbx_0 = model.dim_nbx_0;
    ocp_json.dims.nbu = model.dim_nbu;
    ocp_json.dims.ng = model.dim_ng;
    ocp_json.dims.nh = model.dim_nh;
    if isfield(model, 'dim_ns')
        ocp_json.dims.ns = model.dim_ns;
    end
    if isfield(model, 'dim_nsbx')
        ocp_json.dims.nsbx = model.dim_nsbx;
    end
    if isfield(model, 'dim_nsbu')
        ocp_json.dims.nsbu = model.dim_nsbu;
    end

    % missing in template
    % ocp_json.dims.nsg = model.nsg;

    % missing in MEX
    % ocp_json.dims.npd = model.npd;
    % ocp_json.dims.npd_e = model.npd_e;

    if isfield(model, 'dim_nsh')
        ocp_json.dims.nsh = model.dim_nsh;
    end

    % terminal
    if isfield(model, 'dim_ng_e')
        ocp_json.dims.ng_e = model.dim_ng_e;
    end
    if isfield(model, 'dim_ny_e')
        ocp_json.dims.ny_e = model.dim_ny_e;
    end
    if isfield(model, 'dim_nh_e')
        ocp_json.dims.nh_e = model.dim_nh_e;
    end
    if isfield(model, 'dim_ns_e')
        ocp_json.dims.ns_e = model.dim_ns_e;
    end
    if isfield(model, 'dim_nsh_e')
        ocp_json.dims.nsh_e = model.dim_nsh_e;
    end
    % missing in MEX
    % ocp_json.dims.nbx_e = model.dim_nbx_e;
    % ocp_json.dims.ns_e = model.dim_ns_e;

    %% types
    ocp_json.cost.cost_type = upper(model.cost_type);
    ocp_json.cost.cost_type_e = upper(model.cost_type_e);
    ocp_json.constraints.constr_type = upper(model.constr_type);
    ocp_json.constraints.constr_type_e = upper(model.constr_type_e);

    % parameters
    if model.dim_np > 0
        % TODO: add option to initialize parameters in model.
        warning(['model parameters value cannot be defined (yet) for ocp json.', ...
                    10 'Using zeros(np,1) by default.' 10 'You can update them later using the solver object.']);
        ocp_json.constraints.p = zeros(model.dim_np,1);
    end

    %% constraints
    % initial
    if isfield(model, 'constr_lbx_0')
        ocp_json.constraints.lbx_0 = model.constr_lbx_0;
    else
        warning('missing: constr_lbx_0, using zeros of appropriate dimension.');
        ocp_json.constraints.lbx_0 = zeros(ocp_json.dims.nbx_0, 1);
    end

    if isfield(model, 'constr_ubx_0')
        ocp_json.constraints.ubx_0 = model.constr_ubx_0;
    else
        warning('missing: constr_ubx_0, using zeros of appropriate dimension.');
        ocp_json.constraints.ubx_0 = zeros(ocp_json.dims.nbx_0, 1);
    end
    if isfield(model, 'constr_Jbx_0')
        ocp_json.constraints.idxbx_0 = J_to_idx( model.constr_Jbx_0 );
    end

    % path
    if ocp_json.dims.nbx > 0
        ocp_json.constraints.idxbx = J_to_idx( model.constr_Jbx );
        ocp_json.constraints.lbx = model.constr_lbx;
        ocp_json.constraints.ubx = model.constr_ubx;
    end

    if ocp_json.dims.nbu > 0
        ocp_json.constraints.idxbu = J_to_idx( model.constr_Jbu );
        ocp_json.constraints.lbu = model.constr_lbu;
        ocp_json.constraints.ubu = model.constr_ubu;
    end

    if ocp_json.dims.ng > 0
        ocp_json.constraints.C = model.constr_C;
        ocp_json.constraints.D = model.constr_D;
        ocp_json.constraints.lg = model.constr_lg;
        ocp_json.constraints.ug = model.constr_ug;
    end

    if ocp_json.dims.nh > 0
        ocp_json.con_h.expr = model.constr_expr_h;
        ocp_json.constraints.lh = model.constr_lh;
        ocp_json.constraints.uh = model.constr_uh;
        % TODO(oj): can we get rid of the following?
        if isfield(model, 'sym_x')
            ocp_json.con_h.x = model.sym_x;
        end
        if isfield(model, 'sym_u')
            ocp_json.con_h.u = model.sym_u;
        end
        if isfield(model, 'sym_z')
            ocp_json.con_h.z = model.sym_z;
        end
    end

    if ocp_json.dims.nsbx > 0
        ocp_json.constraints.idxsbx = J_to_idx_slack(model.constr_Jsbx);
        if isfield(model, 'constr_lsbx')
            ocp_json.constraints.lsbx = model.constr_lsbx;
        else
            ocp_json.constraints.lsbx = zeros(ocp_json.dims.nsbx, 1);
        end
        if isfield(model, 'constr_usbx')
            ocp_json.constraints.usbx = model.constr_usbx;
        else
            ocp_json.constraints.usbx = zeros(ocp_json.dims.nsbx, 1);
        end
        % TODO(oj): add nsbx_e properly in Matlab:
        ocp_json.dims.nsbx_e = model.dim_nsbx;
        ocp_json.constraints.idxsbx_e = ocp_json.constraints.idxsbx;
        ocp_json.constraints.lsbx_e = ocp_json.constraints.lsbx;
        ocp_json.constraints.usbx_e = ocp_json.constraints.usbx;
    end


    if ocp_json.dims.nsbu > 0
        ocp_json.constraints.idxsbu = J_to_idx_slack(model.Jsbu);
        if isfield(model, 'constr_lsbu')
            ocp_json.constraints.lsbu = model.constr_lsbu;
        else
            ocp_json.constraints.lsbu = zeros(ocp_json.dims.nsbu, 1);
        end
        if isfield(model, 'constr_usbu')
            ocp_json.constraints.usbu = model.constr_usbu;
        else
            ocp_json.constraints.usbu = zeros(ocp_json.dims.nsbu, 1);
        end
    end

    if ocp_json.dims.nsh > 0
        ocp_json.constraints.idxsh = J_to_idx_slack(model.constr_Jsh);
        if isfield(model, 'constr_lsh')
            ocp_json.constraints.lsh = model.constr_lsh;
        else
            ocp_json.constraints.lsh = zeros(ocp_json.dims.nsh, 1);
        end
        if isfield(model, 'constr_ush')
            ocp_json.constraints.ush = model.constr_ush;
        else
            ocp_json.constraints.ush = zeros(ocp_json.dims.nsh, 1);
        end
    end

    if isfield(model, 'dim_nsg') && model.dim_nsg > 0
        error('dim_nsg > 0 not implmented in code-gen backend');
        % TODO set Jsg
    end


    % terminal
    if ocp_json.dims.ng_e > 0
        ocp_json.constraints.C_e = model.constr_C_e;
        ocp_json.constraints.lg_e = model.constr_lg_e;
        ocp_json.constraints.ug_e = model.constr_ug_e;
    end

    if ocp_json.dims.nh_e > 0    
        ocp_json.con_h_e.expr = model.constr_expr_h_e;
        ocp_json.constraints.lh_e = model.constr_lh_e;
        ocp_json.constraints.uh_e = model.constr_uh_e;
        % TODO(oj): can we get rid of the following?
        if isfield(model, 'sym_x')
            ocp_json.con_h_e.x = model.sym_x;
        end
        if isfield(model, 'sym_u')
            ocp_json.con_h_e.u = model.sym_u;
        end
        if isfield(model, 'sym_z')
            ocp_json.con_h_e.z = model.sym_z;
        end
    end

    if isfield(model, 'dim_nsg_e') && model.dim_nsg_e > 0
        error('dim_nsg_e > 0 not implmented in templating backend');
        % TODO set Jsg_e
    end


    if ocp_json.dims.nsh_e > 0
        ocp_json.constraints.idxsh_e = J_to_idx_slack(model.constr_Jsh_e);
        if isfield(model, 'constr_lsh_e')
            ocp_json.constraints.lsh_e = model.constr_lsh_e;
        else
            ocp_json.constraints.lsh_e = zeros(ocp_json.dims.nsh_e, 1);
        end
        if isfield(model, 'constr_ush_e')
            ocp_json.constraints.ush_e = model.constr_ush_e;
        else
            ocp_json.constraints.ush_e = zeros(ocp_json.dims.nsh_e, 1);
        end
    end

    %% Cost
    if strcmp(model.cost_type, 'linear_ls')
        ocp_json.cost.Vu = model.cost_Vu;
        ocp_json.cost.Vx = model.cost_Vx;
        if isfield(model, 'cost_Vz')
            ocp_json.cost.Vz = model.cost_Vz;
        end
    end

    if strcmp(model.cost_type, 'nonlinear_ls') || strcmp(model.cost_type, 'linear_ls')
        ocp_json.cost.W = model.cost_W;
        if isfield(model, 'cost_y_ref')
            ocp_json.cost.yref = model.cost_y_ref;
        else
			warning(['cost_y_ref not defined for ocp json.' 10 'Using zeros(ny,1) by default.']);
            ocp_json.cost.yref = zeros(model.dim_ny,1);
        end
    end

    if strcmp(model.cost_type_e, 'linear_ls')
        ocp_json.cost.Vx_e = model.cost_Vx_e;
    end

    if strcmp(model.cost_type, 'nonlinear_ls') || strcmp(model.cost_type, 'linear_ls')
        ocp_json.cost.W_e = model.cost_W_e;
        if isfield(model, 'cost_y_ref')
            ocp_json.cost.yref_e = model.cost_y_ref_e;
        else
			warning(['cost_y_ref_e not defined for ocp json.' 10 'Using zeros(ny_e,1) by default.']);
            ocp_json.cost.yref_e = zeros(model.dim_ny_e,1);
        end
    end

    if isfield(model, 'cost_Zl')
        ocp_json.cost.Zl = diag(model.cost_Zl);
    end
    if isfield(model, 'cost_Zu')
        ocp_json.cost.Zu = diag(model.cost_Zu);
    end
    if isfield(model, 'cost_zl')
        ocp_json.cost.zl = model.cost_zl;
    end
    if isfield(model, 'cost_zu')
        ocp_json.cost.zu = model.cost_zu;
    end


    if isfield(model, 'cost_Zl_e')
        ocp_json.cost.Zl_e = diag(model.cost_Zl_e);
    end
    if isfield(model, 'cost_Zu_e')
        ocp_json.cost.Zu_e = diag(model.cost_Zu_e);
    end
    if isfield(model, 'cost_zl_e')
        ocp_json.cost.zl_e = model.cost_zl_e;
    end
    if isfield(model, 'cost_zu_e')
        ocp_json.cost.zu_e = model.cost_zu_e;
    end

    %% dynamics
    if strcmp(obj.opts_struct.sim_method, 'erk')
        ocp_json.model.f_expl_expr = model.dyn_expr_f;
    elseif strcmp(obj.opts_struct.sim_method, 'irk')
        ocp_json.model.f_impl_expr = model.dyn_expr_f;
    else
        error(['integrator ', obj.opts_struct.sim_method, ' not support for templating backend.'])
    end
    %  TODO(oj): add gnsf support;

    ocp_json.model.x = model.sym_x;
    ocp_json.model.u = model.sym_u;
    ocp_json.model.z = model.sym_z;
    ocp_json.model.xdot = model.sym_xdot;
    ocp_json.model.p = model.sym_p;

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
    ncol = size_J(2);
    idx = zeros(ncol,1);
    i_idx = 1;
    for i = 1:nrows
        this_idx = find(J(i,:));
        if length(this_idx) == 1
            idx(i_idx) = this_idx - 1; % strore 0-based index
            i_idx = i_idx + 1;
        elseif length(this_idx) > 1
            error(['J_to_idx: Invalid J matrix. Exiting. Found more than one nonzero in row ' num2str(i)]);
        end
        if J(i,this_idx) ~= 1
            error(['J_to_idx: J matrices can only contain 1s, got J(' num2str(i) ', ' num2str(this_idx) ') = ' num2str(J(i,this_idx)) ]);
        end
    end
end
