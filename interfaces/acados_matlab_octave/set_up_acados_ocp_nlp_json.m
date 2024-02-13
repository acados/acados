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

function ocp_json = set_up_acados_ocp_nlp_json(obj, simulink_opts)

    model = obj.model_struct;
    % create
    ocp_json = acados_template_mex.acados_ocp_nlp_json(simulink_opts);

    % general
    ocp_json.dims.N = obj.opts_struct.param_scheme_N;
    ocp_json.solver_options.tf = model.T;

    ocp_json.code_export_directory = fullfile(pwd, 'c_generated_code');

    if isfield(obj.opts_struct, 'Tsim')
        ocp_json.solver_options.Tsim = obj.opts_struct.Tsim;
    else
        ocp_json.solver_options.Tsim = model.T / obj.opts_struct.param_scheme_N; % for templated integrator
    end

    ocp_json.model.name = model.name;
    ocp_json.name = model.name;

    % modules
    ocp_json.solver_options.qp_solver = upper(obj.opts_struct.qp_solver);
    ocp_json.solver_options.integrator_type = upper(obj.opts_struct.sim_method);
    ocp_json.solver_options.nlp_solver_type = upper(obj.opts_struct.nlp_solver);
    ocp_json.solver_options.collocation_type = upper(obj.opts_struct.collocation_type);

    if strcmp(obj.opts_struct.sim_method, 'irk_gnsf')
        ocp_json.solver_options.integrator_type = 'GNSF';
    end

    N = obj.opts_struct.param_scheme_N;
    % options
    if length(obj.opts_struct.sim_method_num_steps) == N
        ocp_json.solver_options.sim_method_num_steps = obj.opts_struct.sim_method_num_steps;
    else
        ocp_json.solver_options.sim_method_num_steps = obj.opts_struct.sim_method_num_steps * ones(1, N);
    end
    if length(obj.opts_struct.sim_method_num_stages) == N
        ocp_json.solver_options.sim_method_num_stages = obj.opts_struct.sim_method_num_stages;
    else
        ocp_json.solver_options.sim_method_num_stages = obj.opts_struct.sim_method_num_stages * ones(1, N);
    end
    if length(obj.opts_struct.sim_method_jac_reuse) == N
        ocp_json.solver_options.sim_method_jac_reuse = obj.opts_struct.sim_method_jac_reuse;
    else
        ocp_json.solver_options.sim_method_jac_reuse = obj.opts_struct.sim_method_jac_reuse * ones(1, N);
    end

    ocp_json.solver_options.sim_method_newton_iter = obj.opts_struct.sim_method_newton_iter;
    ocp_json.solver_options.nlp_solver_max_iter = obj.opts_struct.nlp_solver_max_iter;
    ocp_json.solver_options.nlp_solver_tol_stat = obj.opts_struct.nlp_solver_tol_stat;
    ocp_json.solver_options.nlp_solver_tol_eq = obj.opts_struct.nlp_solver_tol_eq;
    ocp_json.solver_options.nlp_solver_tol_ineq = obj.opts_struct.nlp_solver_tol_ineq;
    ocp_json.solver_options.nlp_solver_tol_comp = obj.opts_struct.nlp_solver_tol_comp;
    ocp_json.solver_options.nlp_solver_ext_qp_res = obj.opts_struct.nlp_solver_ext_qp_res;
    ocp_json.solver_options.nlp_solver_step_length = obj.opts_struct.nlp_solver_step_length;
    ocp_json.solver_options.globalization = upper(obj.opts_struct.globalization);
    ocp_json.solver_options.alpha_min = obj.opts_struct.alpha_min;
    ocp_json.solver_options.alpha_reduction = obj.opts_struct.alpha_reduction;
    ocp_json.solver_options.line_search_use_sufficient_descent = obj.opts_struct.line_search_use_sufficient_descent;
    ocp_json.solver_options.globalization_use_SOC = obj.opts_struct.globalization_use_SOC;
    ocp_json.solver_options.full_step_dual = obj.opts_struct.full_step_dual;
    ocp_json.solver_options.eps_sufficient_descent = obj.opts_struct.eps_sufficient_descent;
    ocp_json.solver_options.qp_solver_ric_alg = obj.opts_struct.qp_solver_ric_alg;
    ocp_json.solver_options.qp_solver_cond_ric_alg = obj.opts_struct.qp_solver_cond_ric_alg;
    if isfield(obj.opts_struct, 'qp_solver_cond_N')
        ocp_json.solver_options.qp_solver_cond_N = obj.opts_struct.qp_solver_cond_N;
    else
        ocp_json.solver_options.qp_solver_cond_N = obj.opts_struct.param_scheme_N;
    end
    ocp_json.solver_options.qp_solver_iter_max = obj.opts_struct.qp_solver_iter_max;
    if isfield(obj.opts_struct, 'qp_solver_tol_stat')
        ocp_json.solver_options.qp_solver_tol_stat = obj.opts_struct.qp_solver_tol_stat;
    end
    if isfield(obj.opts_struct, 'reg_epsilon')
        ocp_json.solver_options.reg_epsilon = obj.opts_struct.reg_epsilon;
    end
    if isfield(obj.opts_struct, 'qp_solver_tol_eq')
        ocp_json.solver_options.qp_solver_tol_eq = obj.opts_struct.qp_solver_tol_eq;
    end
    if isfield(obj.opts_struct, 'qp_solver_tol_ineq')
        ocp_json.solver_options.qp_solver_tol_ineq = obj.opts_struct.qp_solver_tol_ineq;
    end
    if isfield(obj.opts_struct, 'qp_solver_tol_comp')
        ocp_json.solver_options.qp_solver_tol_comp = obj.opts_struct.qp_solver_tol_comp;
    end
    if isfield(obj.opts_struct, 'qp_solver_warm_start')
        ocp_json.solver_options.qp_solver_warm_start = obj.opts_struct.qp_solver_warm_start;
    end
    if isfield(obj.opts_struct, 'nlp_solver_warm_start_first_qp')
        ocp_json.solver_options.nlp_solver_warm_start_first_qp = obj.opts_struct.nlp_solver_warm_start_first_qp;
    end
    ocp_json.solver_options.levenberg_marquardt = obj.opts_struct.levenberg_marquardt;
    %
    if strcmp(obj.opts_struct.nlp_solver_exact_hessian, 'true')
        ocp_json.solver_options.hessian_approx = 'EXACT';
    else
        ocp_json.solver_options.hessian_approx = 'GAUSS_NEWTON';
    end
    ocp_json.solver_options.regularize_method = upper(obj.opts_struct.regularize_method);

    ocp_json.solver_options.exact_hess_dyn = obj.opts_struct.exact_hess_dyn;
    ocp_json.solver_options.exact_hess_cost = obj.opts_struct.exact_hess_cost;
    ocp_json.solver_options.exact_hess_constr = obj.opts_struct.exact_hess_constr;
    ocp_json.solver_options.fixed_hess = obj.opts_struct.fixed_hess;

    ocp_json.solver_options.ext_fun_compile_flags = obj.opts_struct.ext_fun_compile_flags;

    ocp_json.solver_options.time_steps = obj.opts_struct.time_steps;
    ocp_json.solver_options.print_level = obj.opts_struct.print_level;

    %% dims
    % path
    ocp_json.dims.nx = model.dim_nx;
    ocp_json.dims.nu = model.dim_nu;
    ocp_json.dims.nz = model.dim_nz;
    ocp_json.dims.np = model.dim_np;

    if strcmp(model.cost_type, 'ext_cost')
        ocp_json.dims.ny = 0;
    else
        ocp_json.dims.ny = model.dim_ny;
    end
    ocp_json.dims.nbx = model.dim_nbx;
    ocp_json.dims.nbx_0 = model.dim_nbx_0;
    ocp_json.dims.nbu = model.dim_nbu;
    ocp_json.dims.ng = model.dim_ng;
    ocp_json.dims.nh = model.dim_nh;
    ocp_json.dims.nh_0 = model.dim_nh_0;
    ocp_json.dims.nbxe_0 = model.dim_nbxe_0;
    ocp_json.dims.ns = model.dim_ns;
    ocp_json.dims.nsbx = model.dim_nsbx;
    ocp_json.dims.nsbu = model.dim_nsbu;
    ocp_json.dims.nsg = model.dim_nsg;

    % to match python!
    ocp_json.dims.ns_0 = model.dim_nsbu + model.dim_nsg;

    if isfield(model, 'dim_ny_0')
        ocp_json.dims.ny_0 = model.dim_ny_0;
    elseif strcmp(model.cost_type_0, 'ext_cost')
        ocp_json.dims.ny_0 = 0;
    end
    % missing in MEX
    % ocp_json.dims.nphi;
    % ocp_json.dims.nphi_e;

    if isfield(model, 'dim_nsh')
        ocp_json.dims.nsh = model.dim_nsh;
    end

    % terminal
    ocp_json.dims.nbx_e = model.dim_nbx_e;
    ocp_json.dims.ng_e = model.dim_ng_e;
    if isfield(model, 'dim_ny_e')
        ocp_json.dims.ny_e = model.dim_ny_e;
    elseif strcmp(model.cost_type_e, 'ext_cost')
        ocp_json.dims.ny_e = 0;
    end
    ocp_json.dims.nh_e = model.dim_nh_e;
    ocp_json.dims.ns_e = model.dim_ns_e;
    ocp_json.dims.nsh_e = model.dim_nsh_e;
    ocp_json.dims.nsg_e = model.dim_nsg_e;
    ocp_json.dims.nsbx_e = model.dim_nsbx_e;

    if isfield(model, 'dim_gnsf_nx1')
        ocp_json.dims.gnsf_nx1 = model.dim_gnsf_nx1;
        ocp_json.dims.gnsf_nz1 = model.dim_gnsf_nz1;
        ocp_json.dims.gnsf_nout = model.dim_gnsf_nout;
        ocp_json.dims.gnsf_ny = model.dim_gnsf_ny;
        ocp_json.dims.gnsf_nuhat = model.dim_gnsf_nuhat;
    end

    %% types
    if strcmp(model.cost_type, 'ext_cost')
        ocp_json.cost.cost_type = 'EXTERNAL';
    else
        ocp_json.cost.cost_type = upper(model.cost_type);
    end
    if strcmp(model.cost_type_0, 'ext_cost')
        ocp_json.cost.cost_type_0 = 'EXTERNAL';
    else
        ocp_json.cost.cost_type_0 = upper(model.cost_type_0);
    end
    if strcmp(model.cost_type_e, 'ext_cost')
        ocp_json.cost.cost_type_e = 'EXTERNAL';
    else
        ocp_json.cost.cost_type_e = upper(model.cost_type_e);
    end

    ocp_json.cost.cost_ext_fun_type = model.cost_ext_fun_type;
    if strcmp(model.cost_ext_fun_type, 'generic')
        ocp_json.cost.cost_source_ext_cost = model.cost_source_ext_cost;
        ocp_json.cost.cost_function_ext_cost = model.cost_function_ext_cost;
    end
    ocp_json.cost.cost_ext_fun_type_0 = model.cost_ext_fun_type_0;
    if strcmp(model.cost_ext_fun_type_0, 'generic')
        ocp_json.cost.cost_source_ext_cost_0 = model.cost_source_ext_cost_0;
        ocp_json.cost.cost_function_ext_cost_0 = model.cost_function_ext_cost_0;
    end
    ocp_json.cost.cost_ext_fun_type_e = model.cost_ext_fun_type_e;
    if strcmp(model.cost_ext_fun_type_e, 'generic')
        ocp_json.cost.cost_source_ext_cost_e = model.cost_source_ext_cost_e;
        ocp_json.cost.cost_function_ext_cost_e = model.cost_function_ext_cost_e;
    end

    ocp_json.constraints.constr_type = upper(model.constr_type);
    ocp_json.constraints.constr_type_0 = upper(model.constr_type_0);
    ocp_json.constraints.constr_type_e = upper(model.constr_type_e);

    % parameters
    if model.dim_np > 0
        if isempty(obj.opts_struct.parameter_values)
            warning(['opts_struct.parameter_values are not set.', ...
                        10 'Using zeros(np,1) by default.' 10 'You can update them later using the solver object.']);
            ocp_json.parameter_values = zeros(model.dim_np,1);
        else
            ocp_json.parameter_values = obj.opts_struct.parameter_values(:);
        end
    end

    %% constraints
    % initial
    if isfield(model, 'constr_lbx_0')
        ocp_json.constraints.lbx_0 = model.constr_lbx_0;
    elseif ocp_json.dims.nbx_0 > 0
        warning('missing: constr_lbx_0, using zeros of appropriate dimension.');
        ocp_json.constraints.lbx_0 = zeros(ocp_json.dims.nbx_0, 1);
    end

    if isfield(model, 'constr_ubx_0')
        ocp_json.constraints.ubx_0 = model.constr_ubx_0;
    elseif ocp_json.dims.nbx_0 > 0
        warning('missing: constr_ubx_0, using zeros of appropriate dimension.');
        ocp_json.constraints.ubx_0 = zeros(ocp_json.dims.nbx_0, 1);
    end
    if isfield(model, 'constr_Jbx_0')
        ocp_json.constraints.idxbx_0 = J_to_idx( model.constr_Jbx_0 );
    end
    ocp_json.constraints.idxbxe_0 = model.constr_idxbxe_0;


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
        ocp_json.constraints.lh = model.constr_lh;
        ocp_json.constraints.uh = model.constr_uh;
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
    end


    if ocp_json.dims.nsbu > 0
        ocp_json.constraints.idxsbu = J_to_idx_slack(model.constr_Jsbu);
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

    if ocp_json.dims.nsg > 0
        ocp_json.constraints.idxsg = J_to_idx_slack(model.constr_Jsg);
        if isfield(model, 'constr_lsg')
            ocp_json.constraints.lsg = model.constr_lsg;
        else
            ocp_json.constraints.lsg = zeros(ocp_json.dims.nsg, 1);
        end
        if isfield(model, 'constr_usg')
            ocp_json.constraints.usg = model.constr_usg;
        else
            ocp_json.constraints.usg = zeros(ocp_json.dims.nsg, 1);
        end
    end

    % terminal
    if ocp_json.dims.nbx_e > 0
        ocp_json.constraints.idxbx_e = J_to_idx( model.constr_Jbx_e );
        ocp_json.constraints.lbx_e = model.constr_lbx_e;
        ocp_json.constraints.ubx_e = model.constr_ubx_e;
    end

    if ocp_json.dims.ng_e > 0
        ocp_json.constraints.C_e = model.constr_C_e;
        ocp_json.constraints.lg_e = model.constr_lg_e;
        ocp_json.constraints.ug_e = model.constr_ug_e;
    end

    if ocp_json.dims.nh_e > 0
        ocp_json.constraints.lh_e = model.constr_lh_e;
        ocp_json.constraints.uh_e = model.constr_uh_e;
    end

    if ocp_json.dims.nh_0 > 0
        ocp_json.constraints.lh_0 = model.constr_lh_0;
        ocp_json.constraints.uh_0 = model.constr_uh_0;
    end

    if ocp_json.dims.nsbx_e > 0
        ocp_json.constraints.idxsbx_e = J_to_idx_slack(model.constr_Jsbx_e);
        if isfield(model, 'constr_lsbx_e')
            ocp_json.constraints.lsbx_e = model.constr_lsbx_e;
        else
            ocp_json.constraints.lsbx_e = zeros(ocp_json.dims.nsbx_e, 1);
        end
        if isfield(model, 'constr_usbx_e')
            ocp_json.constraints.usbx_e = model.constr_usbx_e;
        else
            ocp_json.constraints.usbx_e = zeros(ocp_json.dims.nsbx_e, 1);
        end
    end

    if ocp_json.dims.nsg_e > 0
        ocp_json.constraints.idxsg_e = J_to_idx_slack(model.constr_Jsg_e);
        if isfield(model, 'constr_lsg_e')
            ocp_json.constraints.lsg_e = model.constr_lsg_e;
        else
            ocp_json.constraints.lsg_e = zeros(ocp_json.dims.nsg_e, 1);
        end
        if isfield(model, 'constr_usg_e')
            ocp_json.constraints.usg_e = model.constr_usg_e;
        else
            ocp_json.constraints.usg_e = zeros(ocp_json.dims.nsg_e, 1);
        end
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

    if strcmp(model.cost_type_0, 'linear_ls')
        ocp_json.cost.Vu_0 = model.cost_Vu_0;
        ocp_json.cost.Vx_0 = model.cost_Vx_0;
        if isfield(model, 'cost_Vz_0')
            ocp_json.cost.Vz_0 = model.cost_Vz_0;
        end
    end
    if strcmp(model.cost_type_0, 'nonlinear_ls') || strcmp(model.cost_type_0, 'linear_ls')
        ocp_json.cost.W_0 = model.cost_W_0;
        if isfield(model, 'cost_y_ref_0')
            ocp_json.cost.yref_0 = model.cost_y_ref_0;
        else
			warning(['cost_y_ref_0 not defined for ocp json.' 10 'Using zeros(ny_0,1) by default.']);
            ocp_json.cost.yref_0 = zeros(model.dim_ny_0,1);
        end
    end

    if isfield(model, 'cost_Vx_e')
        ocp_json.cost.Vx_e = model.cost_Vx_e;
    end

    if strcmp(model.cost_type_e, 'nonlinear_ls') || strcmp(model.cost_type_e, 'linear_ls')
        if isfield(model, 'cost_W_e')
            ocp_json.cost.W_e = model.cost_W_e;
        end
        if isfield(model, 'cost_y_ref_e')
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

    if isfield(model, 'cost_Zl_0')
        ocp_json.cost.Zl_0 = diag(model.cost_Zl_0);
    end
    if isfield(model, 'cost_Zu_0')
        ocp_json.cost.Zu_0 = diag(model.cost_Zu_0);
    end
    if isfield(model, 'cost_zl_0')
        ocp_json.cost.zl_0 = model.cost_zl_0;
    end
    if isfield(model, 'cost_zu_0')
        ocp_json.cost.zu_0 = model.cost_zu_0;
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
        if strcmp(model.dyn_ext_fun_type, 'casadi')
            ocp_json.model.f_impl_expr = model.dyn_expr_f;
        elseif strcmp(model.dyn_ext_fun_type, 'generic')
            ocp_json.model.dyn_generic_source = model.dyn_generic_source;
        end
    elseif strcmp(obj.opts_struct.sim_method, 'discrete')
        ocp_json.model.dyn_ext_fun_type = model.dyn_ext_fun_type;
        if strcmp(model.dyn_ext_fun_type, 'casadi')
            ocp_json.model.f_phi_expr = model.dyn_expr_phi;
        elseif strcmp(model.dyn_ext_fun_type, 'generic')
            ocp_json.model.dyn_generic_source = model.dyn_generic_source;
            if isfield(model, 'dyn_disc_fun_jac_hess')
                ocp_json.model.dyn_disc_fun_jac_hess = model.dyn_disc_fun_jac_hess;
            end
            if isfield(model, 'dyn_disc_fun_jac')
                ocp_json.model.dyn_disc_fun_jac = model.dyn_disc_fun_jac;
            end
            ocp_json.model.dyn_disc_fun = model.dyn_disc_fun;
        end
    elseif strcmp(obj.opts_struct.sim_method, 'irk_gnsf')
        ocp_json.model.gnsf.A = model.dyn_gnsf_A;
        ocp_json.model.gnsf.B = model.dyn_gnsf_B;
        ocp_json.model.gnsf.C = model.dyn_gnsf_C;
        ocp_json.model.gnsf.E = model.dyn_gnsf_E;
        ocp_json.model.gnsf.c = model.dyn_gnsf_c;
        ocp_json.model.gnsf.A_LO = model.dyn_gnsf_A_LO;
        ocp_json.model.gnsf.B_LO = model.dyn_gnsf_B_LO;
        ocp_json.model.gnsf.E_LO = model.dyn_gnsf_E_LO;
        ocp_json.model.gnsf.c_LO = model.dyn_gnsf_c_LO;

        ocp_json.model.gnsf.L_x = model.dyn_gnsf_L_x;
        ocp_json.model.gnsf.L_u = model.dyn_gnsf_L_u;
        ocp_json.model.gnsf.L_xdot = model.dyn_gnsf_L_xdot;
        ocp_json.model.gnsf.L_z = model.dyn_gnsf_L_z;

        ocp_json.model.gnsf.expr_phi = model.dyn_gnsf_expr_phi;
        ocp_json.model.gnsf.expr_f_lo = model.dyn_gnsf_expr_f_lo;

        ocp_json.model.gnsf.ipiv_x = model.dyn_gnsf_ipiv_x;
        ocp_json.model.gnsf.idx_perm_x = model.dyn_gnsf_idx_perm_x;
        ocp_json.model.gnsf.ipiv_z = model.dyn_gnsf_ipiv_z;
        ocp_json.model.gnsf.idx_perm_z = model.dyn_gnsf_idx_perm_z;
        ocp_json.model.gnsf.ipiv_f = model.dyn_gnsf_ipiv_f;
        ocp_json.model.gnsf.idx_perm_f = model.dyn_gnsf_idx_perm_f;

        ocp_json.model.gnsf.nontrivial_f_LO = model.dyn_gnsf_nontrivial_f_LO;
        ocp_json.model.gnsf.purely_linear = model.dyn_gnsf_purely_linear;

        ocp_json.model.gnsf.y = model.sym_gnsf_y;
        ocp_json.model.gnsf.uhat = model.sym_gnsf_uhat;
    else
        error(['integrator ', obj.opts_struct.sim_method, ' not support for templating backend.'])
    end

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
