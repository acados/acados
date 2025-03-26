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

function ocp = setup_AcadosOcp_from_legacy_ocp_description(model_old, opts_old, simulink_opts)

    model = model_old.model_struct;
    opts_struct = opts_old.opts_struct;

    if nargin < 3
        simulink_opts = get_acados_simulink_opts();
    end

    % create
    ocp = AcadosOcp();
    ocp.simulink_opts = simulink_opts;

    % general
    ocp.solver_options.N_horizon = opts_struct.param_scheme_N;
    ocp.solver_options.tf = model.T;

    ocp.model.name = model.name;

    % modules
    ocp.solver_options.qp_solver = upper(opts_struct.qp_solver);
    ocp.solver_options.nlp_solver_type = upper(opts_struct.nlp_solver);
    ocp.solver_options.collocation_type = upper(opts_struct.collocation_type);

    if strcmp(opts_struct.sim_method, 'irk_gnsf')
        ocp.solver_options.integrator_type = 'GNSF';
    else
        ocp.solver_options.integrator_type = upper(opts_struct.sim_method);
    end

    % options
    ocp.solver_options.sim_method_num_steps = opts_struct.sim_method_num_steps;
    ocp.solver_options.sim_method_num_stages = opts_struct.sim_method_num_stages;
    ocp.solver_options.sim_method_jac_reuse = opts_struct.sim_method_jac_reuse;

    ocp.solver_options.sim_method_newton_iter = opts_struct.sim_method_newton_iter;
    ocp.solver_options.sim_method_newton_tol = opts_struct.sim_method_newton_tol;
    ocp.solver_options.nlp_solver_max_iter = opts_struct.nlp_solver_max_iter;
    ocp.solver_options.nlp_solver_tol_stat = opts_struct.nlp_solver_tol_stat;
    ocp.solver_options.nlp_solver_tol_eq = opts_struct.nlp_solver_tol_eq;
    ocp.solver_options.nlp_solver_tol_ineq = opts_struct.nlp_solver_tol_ineq;
    ocp.solver_options.nlp_solver_tol_comp = opts_struct.nlp_solver_tol_comp;
    ocp.solver_options.nlp_solver_ext_qp_res = opts_struct.nlp_solver_ext_qp_res;
    ocp.solver_options.globalization = upper(opts_struct.globalization);
    ocp.solver_options.globalization_fixed_step_length = opts_struct.globalization_fixed_step_length;
    ocp.solver_options.globalization_alpha_min = opts_struct.globalization_alpha_min;
    ocp.solver_options.globalization_alpha_reduction = opts_struct.globalization_alpha_reduction;
    ocp.solver_options.globalization_line_search_use_sufficient_descent = opts_struct.globalization_line_search_use_sufficient_descent;
    ocp.solver_options.globalization_use_SOC = opts_struct.globalization_use_SOC;
    ocp.solver_options.globalization_full_step_dual = opts_struct.globalization_full_step_dual;
    ocp.solver_options.globalization_eps_sufficient_descent = opts_struct.globalization_eps_sufficient_descent;
    ocp.solver_options.qp_solver_ric_alg = opts_struct.qp_solver_ric_alg;
    ocp.solver_options.qp_solver_cond_ric_alg = opts_struct.qp_solver_cond_ric_alg;
    ocp.solver_options.qp_solver_mu0 = opts_struct.qp_solver_mu0;
    ocp.solver_options.store_iterates = opts_struct.store_iterates;
    ocp.json_file = opts_struct.json_file;
    if isfield(opts_struct, 'qp_solver_cond_N')
        ocp.solver_options.qp_solver_cond_N = opts_struct.qp_solver_cond_N;
    else
        ocp.solver_options.qp_solver_cond_N = opts_struct.param_scheme_N;
    end

    ocp.solver_options.qp_solver_iter_max = opts_struct.qp_solver_iter_max;
    if isfield(opts_struct, 'qp_solver_tol_stat')
        ocp.solver_options.qp_solver_tol_stat = opts_struct.qp_solver_tol_stat;
    end
    if isfield(opts_struct, 'reg_epsilon')
        ocp.solver_options.reg_epsilon = opts_struct.reg_epsilon;
    end
    if isfield(opts_struct, 'qp_solver_tol_eq')
        ocp.solver_options.qp_solver_tol_eq = opts_struct.qp_solver_tol_eq;
    end
    if isfield(opts_struct, 'qp_solver_tol_ineq')
        ocp.solver_options.qp_solver_tol_ineq = opts_struct.qp_solver_tol_ineq;
    end
    if isfield(opts_struct, 'qp_solver_tol_comp')
        ocp.solver_options.qp_solver_tol_comp = opts_struct.qp_solver_tol_comp;
    end
    if isfield(opts_struct, 'qp_solver_warm_start')
        ocp.solver_options.qp_solver_warm_start = opts_struct.qp_solver_warm_start;
    end
    if isfield(opts_struct, 'nlp_solver_warm_start_first_qp')
        ocp.solver_options.nlp_solver_warm_start_first_qp = opts_struct.nlp_solver_warm_start_first_qp;
    end
    ocp.solver_options.levenberg_marquardt = opts_struct.levenberg_marquardt;
    %
    if strcmp(opts_struct.nlp_solver_exact_hessian, 'true')
        ocp.solver_options.hessian_approx = 'EXACT';
    else
        ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
    end
    ocp.solver_options.regularize_method = upper(opts_struct.regularize_method);

    ocp.solver_options.exact_hess_dyn = opts_struct.exact_hess_dyn;
    ocp.solver_options.exact_hess_cost = opts_struct.exact_hess_cost;
    ocp.solver_options.exact_hess_constr = opts_struct.exact_hess_constr;
    ocp.solver_options.fixed_hess = opts_struct.fixed_hess;

    ocp.solver_options.ext_fun_compile_flags = opts_struct.ext_fun_compile_flags;
    ocp.solver_options.ext_fun_expand_dyn = opts_struct.ext_fun_expand_dyn;
    ocp.solver_options.ext_fun_expand_cost = opts_struct.ext_fun_expand_cost;
    ocp.solver_options.ext_fun_expand_constr = opts_struct.ext_fun_expand_constr;
    ocp.solver_options.ext_fun_expand_precompute = opts_struct.ext_fun_expand_precompute;

    ocp.solver_options.time_steps = opts_struct.time_steps;
    ocp.solver_options.shooting_nodes = opts_struct.shooting_nodes;
    ocp.solver_options.print_level = opts_struct.print_level;

    ocp.solver_options.timeout_max_time = opts_struct.timeout_max_time;
    ocp.solver_options.timeout_heuristic = opts_struct.timeout_heuristic;

    % compile mex interface (without model dependency)
    if strcmp(opts_struct.compile_interface, 'true')
        ocp.solver_options.compile_interface = true;
    elseif strcmp(opts_struct.compile_interface, 'false')
        ocp.solver_options.compile_interface = false;
    elseif strcmp(opts_struct.compile_interface, 'auto')
        ocp.solver_options.compile_interface = [];
    else
        error(['Unknown value for `compile_interface`' opts_struct.compile_interface]);
    end

    %% types
    if strcmp(model.cost_type, 'ext_cost')
        ocp.cost.cost_type = 'EXTERNAL';
    else
        ocp.cost.cost_type = upper(model.cost_type);
    end
    if strcmp(model.cost_type_0, 'ext_cost')
        ocp.cost.cost_type_0 = 'EXTERNAL';
    else
        ocp.cost.cost_type_0 = upper(model.cost_type_0);
    end
    if strcmp(model.cost_type_e, 'ext_cost')
        ocp.cost.cost_type_e = 'EXTERNAL';
    else
        ocp.cost.cost_type_e = upper(model.cost_type_e);
    end

    ocp.constraints.constr_type = upper(model.constr_type);
    ocp.constraints.constr_type_0 = upper(model.constr_type_0);
    ocp.constraints.constr_type_e = upper(model.constr_type_e);

    % parameters
    ocp.parameter_values = opts_struct.parameter_values;
    ocp.p_global_values = opts_struct.p_global_values;

    %% constraints
    constraints_fields_map = struct(...
        'constr_idxbxe_0', 'idxbxe_0', ...
        'constr_lbx_0', 'lbx_0', ...
        'constr_ubx_0', 'ubx_0', ...
        'constr_Jbx_0', 'idxbx_0', ...
        'constr_lh_0', 'lh_0', ...
        'constr_uh_0', 'uh_0', ...
        'constr_lsh_0', 'lsh_0', ...
        'constr_ush_0', 'ush_0', ...
        'constr_Jsh_0', 'idxsh_0', ...
        'constr_lbu', 'lbu', ...
        'constr_ubu', 'ubu', ...
        'constr_Jbu', 'idxbu', ...
        'constr_lbx', 'lbx', ...
        'constr_ubx', 'ubx', ...
        'constr_Jbx', 'idxbx', ...
        'constr_C', 'C', ...
        'constr_D', 'D', ...
        'constr_lg', 'lg', ...
        'constr_ug', 'ug', ...
        'constr_lh', 'lh', ...
        'constr_uh', 'uh', ...
        'constr_lsbx', 'lsbx', ...
        'constr_usbx', 'usbx', ...
        'constr_Jsbx', 'idxsbx', ...
        'constr_lsbu', 'lsbu', ...
        'constr_usbu', 'usbu', ...
        'constr_Jsbu', 'idxsbu', ...
        'constr_lsh', 'lsh', ...
        'constr_ush', 'ush', ...
        'constr_Jsh', 'idxsh', ...
        'constr_lsg', 'lsg', ...
        'constr_usg', 'usg', ...
        'constr_Jsg', 'idxsg', ...
        'constr_lbx_e', 'lbx_e', ...
        'constr_ubx_e', 'ubx_e', ...
        'constr_Jbx_e', 'idxbx_e', ...
        'constr_C_e', 'C_e', ...
        'constr_lg_e', 'lg_e', ...
        'constr_ug_e', 'ug_e', ...
        'constr_lh_e', 'lh_e', ...
        'constr_uh_e', 'uh_e', ...
        'constr_lsbx_e', 'lsbx_e', ...
        'constr_usbx_e', 'usbx_e', ...
        'constr_Jsbx_e', 'idxsbx_e', ...
        'constr_lsg_e', 'lsg_e', ...
        'constr_usg_e', 'usg_e', ...
        'constr_Jsg_e', 'idxsg_e', ...
        'constr_lsh_e', 'lsh_e', ...
        'constr_ush_e', 'ush_e', ...
        'constr_Jsh_e', 'idxsh_e' ...
        );
    fields = fieldnames(constraints_fields_map);
    num_fields = length(fields);
    for n=1:num_fields
        old_field = fields{n};
        if isfield(model, old_field)
        new_field = constraints_fields_map.(old_field);

            if strncmp(old_field, 'constr_Jb', 9)
                ocp.constraints.(new_field) = J_to_idx(model.(old_field));
            elseif strncmp(old_field, 'constr_Js', 9)
                ocp.constraints.(new_field) = J_to_idx_slack(model.(old_field));
            else
                ocp.constraints.(new_field) = model.(old_field);
            end
        end
    end

    model_fields_map = struct( ...
        'constr_expr_h_0', 'con_h_expr_0', ...
        'constr_expr_h', 'con_h_expr', ...
        'constr_expr_h_e', 'con_h_expr_e', ...
        'sym_x', 'x', ...
        'sym_xdot', 'xdot', ...
        'sym_u', 'u', ...
        'sym_p', 'p', ...
        'sym_p_global', 'p_global', ...
        'sym_z', 'z', ...
        'sym_t', 't', ...
        'cost_expr_ext_cost_0', 'cost_expr_ext_cost_0', ...
        'cost_expr_ext_cost', 'cost_expr_ext_cost', ...
        'cost_expr_ext_cost_e', 'cost_expr_ext_cost_e', ...
        'cost_expr_y_0', 'cost_y_expr_0', ...
        'cost_expr_y', 'cost_y_expr', ...
        'cost_expr_y_e', 'cost_y_expr_e', ...
        'dyn_ext_fun_type', 'dyn_ext_fun_type' ...
        );

    fields = fieldnames(model_fields_map);
    num_fields = length(fields);
    for n=1:num_fields
        old_field = fields{n};
        if isfield(model, old_field)
            new_field = model_fields_map.(old_field);
            ocp.model.(new_field) = model.(old_field);
        end
    end

    %% Cost
    cost_fields_map = struct( ...
        'cost_ext_fun_type_0', 'cost_ext_fun_type_0', ...
        'cost_ext_fun_type', 'cost_ext_fun_type', ...
        'cost_ext_fun_type_e', 'cost_ext_fun_type_e', ...
        'cost_source_ext_cost_0', 'cost_source_ext_cost_0', ...
        'cost_source_ext_cost', 'cost_source_ext_cost', ...
        'cost_source_ext_cost_e', 'cost_source_ext_cost_e', ...
        'cost_function_ext_cost_0', 'cost_function_ext_cost_0', ...
        'cost_function_ext_cost', 'cost_function_ext_cost', ...
        'cost_function_ext_cost_e', 'cost_function_ext_cost_e', ...
        'cost_W_0', 'W_0', ...
        'cost_W', 'W', ...
        'cost_W_e', 'W_e', ...
        'cost_y_ref_0', 'yref_0', ...
        'cost_y_ref', 'yref', ...
        'cost_y_ref_e', 'yref_e', ...
        'cost_Vx_0', 'Vx_0', ...
        'cost_Vx', 'Vx', ...
        'cost_Vx_e', 'Vx_e', ...
        'cost_Vu_0', 'Vu_0', ...
        'cost_Vu', 'Vu', ...
        'cost_Vz_0', 'Vz_0', ...
        'cost_Vz', 'Vz', ...
        'cost_Zl', 'Zl', ...
        'cost_Zu', 'Zu', ...
        'cost_zl', 'zl', ...
        'cost_zu', 'zu', ...
        'cost_Zl_0', 'Zl_0', ...
        'cost_Zu_0', 'Zu_0', ...
        'cost_zl_0', 'zl_0', ...
        'cost_zu_0', 'zu_0', ...
        'cost_Zl_e', 'Zl_e', ...
        'cost_Zu_e', 'Zu_e', ...
        'cost_zl_e', 'zl_e', ...
        'cost_zu_e', 'zu_e' ...
        );
    fields = fieldnames(cost_fields_map);
    num_fields = length(fields);
    for n=1:num_fields
        old_field = fields{n};
        if isfield(model, old_field)
            new_field = cost_fields_map.(old_field);
            ocp.cost.(new_field) = model.(old_field);
            if strncmp(new_field, 'Z', 1)
                ocp.cost.(new_field) = diag(model.(old_field));
            end
        end
    end

    %% dynamics
    if strcmp(opts_struct.sim_method, 'erk')
        ocp.model.f_expl_expr = model.dyn_expr_f;
    elseif strcmp(opts_struct.sim_method, 'irk') || strcmp(opts_struct.sim_method, 'irk_gnsf')
        if strcmp(model.dyn_ext_fun_type, 'casadi')
            ocp.model.f_impl_expr = model.dyn_expr_f;
        elseif strcmp(model.dyn_ext_fun_type, 'generic')
            ocp.model.dyn_generic_source = model.dyn_generic_source;
        end
    elseif strcmp(opts_struct.sim_method, 'discrete')
        ocp.model.dyn_ext_fun_type = model.dyn_ext_fun_type;
        if strcmp(model.dyn_ext_fun_type, 'casadi')
            ocp.model.disc_dyn_expr = model.dyn_expr_phi;
        elseif strcmp(model.dyn_ext_fun_type, 'generic')
            ocp.model.dyn_generic_source = model.dyn_generic_source;
            if isfield(model, 'dyn_disc_fun_jac_hess')
                ocp.model.dyn_disc_fun_jac_hess = model.dyn_disc_fun_jac_hess;
            end
            if isfield(model, 'dyn_disc_fun_jac')
                ocp.model.dyn_disc_fun_jac = model.dyn_disc_fun_jac;
            end
            ocp.model.dyn_disc_fun = model.dyn_disc_fun;
        end
    elseif strcmp(opts_struct.sim_method, 'irk_gnsf')
        gnsf_fields_map = struct(...
        'dyn_gnsf_A', 'dyn_gnsf_A', ...
        'dyn_gnsf_B', 'dyn_gnsf_B', ...
        'dyn_gnsf_C', 'dyn_gnsf_C', ...
        'dyn_gnsf_E', 'dyn_gnsf_E', ...
        'dyn_gnsf_c', 'dyn_gnsf_c', ...
        'dyn_gnsf_A_LO', 'dyn_gnsf_A_LO', ...
        'dyn_gnsf_B_LO', 'dyn_gnsf_B_LO', ...
        'dyn_gnsf_E_LO', 'dyn_gnsf_E_LO', ...
        'dyn_gnsf_c_LO', 'dyn_gnsf_c_LO', ...
        'dyn_gnsf_L_x', 'dyn_gnsf_L_x', ...
        'dyn_gnsf_L_u', 'dyn_gnsf_L_u', ...
        'dyn_gnsf_L_xdot', 'dyn_gnsf_L_xdot', ...
        'dyn_gnsf_L_z', 'dyn_gnsf_L_z', ...
        'dyn_gnsf_expr_phi', 'dyn_gnsf_expr_phi', ...
        'dyn_gnsf_expr_f_lo', 'dyn_gnsf_expr_f_lo', ...
        'dyn_gnsf_ipiv_x', 'dyn_gnsf_ipiv_x', ...
        'dyn_gnsf_idx_perm_x', 'dyn_gnsf_idx_perm_x', ...
        'dyn_gnsf_ipiv_z', 'dyn_gnsf_ipiv_z', ...
        'dyn_gnsf_idx_perm_z', 'dyn_gnsf_idx_perm_z', ...
        'dyn_gnsf_ipiv_f', 'dyn_gnsf_ipiv_f', ...
        'dyn_gnsf_idx_perm_f', 'dyn_gnsf_idx_perm_f', ...
        'dyn_gnsf_nontrivial_f_LO', 'dyn_gnsf_nontrivial_f_LO', ...
        'dyn_gnsf_purely_linear', 'dyn_gnsf_purely_linear', ...
        'sym_gnsf_y', 'sym_gnsf_y', ...
        'sym_gnsf_uhat', 'sym_gnsf_uhat' ...
        );
        fields = fieldnames(gnsf_fields_map);
        num_fields = length(fields);
        for n=1:num_fields
            old_field = fields{n};
            if isfield(model, old_field)
                new_field = gnsf_fields_map.(old_field);
                ocp.model.(new_field) = model.(old_field);
            end
        end
    else
        error(['integrator ', opts_struct.sim_method, ' not support for templating backend.'])
    end

    % check if has_x0 is true. If true inequalities are set to equalities.
    % NOTE: This adds an additional call to make_consistent, but is
    % necessary for the check for has_x0 as we need nx and nbx_0 available,
    % as well as consistent dimensions for lbx_0 and ubx_0.
    ocp.model.make_consistent(ocp.dims);
    if ocp.dims.nx == ocp.dims.nbx_0 && all(sort(reshape(ocp.constraints.idxbx_0, 1, [])) == 0:(ocp.dims.nx-1)) && all(ocp.constraints.lbx_0 == ocp.constraints.ubx_0)
        ocp.constraints.x0 = ocp.constraints.lbx_0; % this will set has_x0 to true and set additional fields
    end
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
            idx(i_idx) = i - 1; % store 0-based index
            i_idx = i_idx + 1;
        elseif length(this_idx) > 1
            error(['J_to_idx_slack: Invalid J matrix. Exiting. Found more than one nonzero in row ' num2str(i)]);
        end
        if J(i,this_idx) ~= 1
            error(['J_to_idx_slack: J matrices can only contain 1s, got J(' num2str(i) ', ' num2str(this_idx) ') = ' num2str(J(i,this_idx)) ]);
        end
    end
    if i_idx ~= ncol + 1
        error('J_to_idx_slack: J must contain a 1 in every column!')
    end
end
