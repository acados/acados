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

function ocp = setup_ocp(obj, simulink_opts)

    model = obj.model_struct;
    % create
    ocp = acados_template_mex.AcadosOcp(simulink_opts);

    % general
    ocp.dims.N = obj.opts_struct.param_scheme_N;
    ocp.solver_options.tf = model.T;

    ocp.code_export_directory = fullfile(pwd, 'c_generated_code');

    if isfield(obj.opts_struct, 'Tsim')
        ocp.solver_options.Tsim = obj.opts_struct.Tsim;
    else
        ocp.solver_options.Tsim = model.T / obj.opts_struct.param_scheme_N; % for templated integrator
    end

    ocp.model.name = model.name;
    ocp.name = model.name;

    % modules
    ocp.solver_options.qp_solver = upper(obj.opts_struct.qp_solver);
    ocp.solver_options.integrator_type = upper(obj.opts_struct.sim_method);
    ocp.solver_options.nlp_solver_type = upper(obj.opts_struct.nlp_solver);
    ocp.solver_options.collocation_type = upper(obj.opts_struct.collocation_type);

    if strcmp(obj.opts_struct.sim_method, 'irk_gnsf')
        ocp.solver_options.integrator_type = 'GNSF';
    end

    N = obj.opts_struct.param_scheme_N;
    % options
    if length(obj.opts_struct.sim_method_num_steps) == N
        ocp.solver_options.sim_method_num_steps = obj.opts_struct.sim_method_num_steps;
    else
        ocp.solver_options.sim_method_num_steps = obj.opts_struct.sim_method_num_steps * ones(1, N);
    end
    if length(obj.opts_struct.sim_method_num_stages) == N
        ocp.solver_options.sim_method_num_stages = obj.opts_struct.sim_method_num_stages;
    else
        ocp.solver_options.sim_method_num_stages = obj.opts_struct.sim_method_num_stages * ones(1, N);
    end
    if length(obj.opts_struct.sim_method_jac_reuse) == N
        ocp.solver_options.sim_method_jac_reuse = obj.opts_struct.sim_method_jac_reuse;
    else
        ocp.solver_options.sim_method_jac_reuse = obj.opts_struct.sim_method_jac_reuse * ones(1, N);
    end

    ocp.solver_options.sim_method_newton_iter = obj.opts_struct.sim_method_newton_iter;
    ocp.solver_options.nlp_solver_max_iter = obj.opts_struct.nlp_solver_max_iter;
    ocp.solver_options.nlp_solver_tol_stat = obj.opts_struct.nlp_solver_tol_stat;
    ocp.solver_options.nlp_solver_tol_eq = obj.opts_struct.nlp_solver_tol_eq;
    ocp.solver_options.nlp_solver_tol_ineq = obj.opts_struct.nlp_solver_tol_ineq;
    ocp.solver_options.nlp_solver_tol_comp = obj.opts_struct.nlp_solver_tol_comp;
    ocp.solver_options.nlp_solver_ext_qp_res = obj.opts_struct.nlp_solver_ext_qp_res;
    ocp.solver_options.nlp_solver_step_length = obj.opts_struct.nlp_solver_step_length;
    ocp.solver_options.globalization = upper(obj.opts_struct.globalization);
    ocp.solver_options.alpha_min = obj.opts_struct.alpha_min;
    ocp.solver_options.alpha_reduction = obj.opts_struct.alpha_reduction;
    ocp.solver_options.line_search_use_sufficient_descent = obj.opts_struct.line_search_use_sufficient_descent;
    ocp.solver_options.globalization_use_SOC = obj.opts_struct.globalization_use_SOC;
    ocp.solver_options.full_step_dual = obj.opts_struct.full_step_dual;
    ocp.solver_options.eps_sufficient_descent = obj.opts_struct.eps_sufficient_descent;
    ocp.solver_options.qp_solver_ric_alg = obj.opts_struct.qp_solver_ric_alg;
    ocp.solver_options.qp_solver_cond_ric_alg = obj.opts_struct.qp_solver_cond_ric_alg;
    if isfield(obj.opts_struct, 'qp_solver_cond_N')
        ocp.solver_options.qp_solver_cond_N = obj.opts_struct.qp_solver_cond_N;
    else
        ocp.solver_options.qp_solver_cond_N = obj.opts_struct.param_scheme_N;
    end
    ocp.solver_options.qp_solver_iter_max = obj.opts_struct.qp_solver_iter_max;
    if isfield(obj.opts_struct, 'qp_solver_tol_stat')
        ocp.solver_options.qp_solver_tol_stat = obj.opts_struct.qp_solver_tol_stat;
    end
    if isfield(obj.opts_struct, 'reg_epsilon')
        ocp.solver_options.reg_epsilon = obj.opts_struct.reg_epsilon;
    end
    if isfield(obj.opts_struct, 'qp_solver_tol_eq')
        ocp.solver_options.qp_solver_tol_eq = obj.opts_struct.qp_solver_tol_eq;
    end
    if isfield(obj.opts_struct, 'qp_solver_tol_ineq')
        ocp.solver_options.qp_solver_tol_ineq = obj.opts_struct.qp_solver_tol_ineq;
    end
    if isfield(obj.opts_struct, 'qp_solver_tol_comp')
        ocp.solver_options.qp_solver_tol_comp = obj.opts_struct.qp_solver_tol_comp;
    end
    if isfield(obj.opts_struct, 'qp_solver_warm_start')
        ocp.solver_options.qp_solver_warm_start = obj.opts_struct.qp_solver_warm_start;
    end
    if isfield(obj.opts_struct, 'nlp_solver_warm_start_first_qp')
        ocp.solver_options.nlp_solver_warm_start_first_qp = obj.opts_struct.nlp_solver_warm_start_first_qp;
    end
    ocp.solver_options.levenberg_marquardt = obj.opts_struct.levenberg_marquardt;
    %
    if strcmp(obj.opts_struct.nlp_solver_exact_hessian, 'true')
        ocp.solver_options.hessian_approx = 'EXACT';
    else
        ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
    end
    ocp.solver_options.regularize_method = upper(obj.opts_struct.regularize_method);

    ocp.solver_options.exact_hess_dyn = obj.opts_struct.exact_hess_dyn;
    ocp.solver_options.exact_hess_cost = obj.opts_struct.exact_hess_cost;
    ocp.solver_options.exact_hess_constr = obj.opts_struct.exact_hess_constr;
    ocp.solver_options.fixed_hess = obj.opts_struct.fixed_hess;

    ocp.solver_options.ext_fun_compile_flags = obj.opts_struct.ext_fun_compile_flags;

    ocp.solver_options.time_steps = obj.opts_struct.time_steps;
    ocp.solver_options.print_level = obj.opts_struct.print_level;

    %% dims
    % path
    ocp.dims.nx = model.dim_nx;
    ocp.dims.nu = model.dim_nu;
    ocp.dims.nz = model.dim_nz;
    ocp.dims.np = model.dim_np;

    if strcmp(model.cost_type, 'ext_cost')
        ocp.dims.ny = 0;
    else
        ocp.dims.ny = model.dim_ny;
    end
    ocp.dims.nbx = model.dim_nbx;
    ocp.dims.nbx_0 = model.dim_nbx_0;
    ocp.dims.nbu = model.dim_nbu;
    ocp.dims.ng = model.dim_ng;
    ocp.dims.nh = model.dim_nh;
    ocp.dims.nh_0 = model.dim_nh_0;
    ocp.dims.nbxe_0 = model.dim_nbxe_0;
    ocp.dims.ns = model.dim_ns;
    ocp.dims.nsbx = model.dim_nsbx;
    ocp.dims.nsbu = model.dim_nsbu;
    ocp.dims.nsg = model.dim_nsg;
    ocp.dims.ns_0 = model.dim_ns_0;
    ocp.dims.nsh_0 = model.dim_nsh_0;
    ocp.dims.nsphi_0 = model.dim_nsphi_0;

    if isfield(model, 'dim_ny_0')
        ocp.dims.ny_0 = model.dim_ny_0;
    elseif strcmp(model.cost_type_0, 'ext_cost')
        ocp.dims.ny_0 = 0;
    end
    % missing in MEX
    % ocp.dims.nphi;
    % ocp.dims.nphi_e;

    if isfield(model, 'dim_nsh')
        ocp.dims.nsh = model.dim_nsh;
    end

    % terminal
    ocp.dims.nbx_e = model.dim_nbx_e;
    ocp.dims.ng_e = model.dim_ng_e;
    if isfield(model, 'dim_ny_e')
        ocp.dims.ny_e = model.dim_ny_e;
    elseif strcmp(model.cost_type_e, 'ext_cost')
        ocp.dims.ny_e = 0;
    end
    ocp.dims.nh_e = model.dim_nh_e;
    ocp.dims.ns_e = model.dim_ns_e;
    ocp.dims.nsh_e = model.dim_nsh_e;
    ocp.dims.nsg_e = model.dim_nsg_e;
    ocp.dims.nsbx_e = model.dim_nsbx_e;

    if isfield(model, 'dim_gnsf_nx1')
        ocp.dims.gnsf_nx1 = model.dim_gnsf_nx1;
        ocp.dims.gnsf_nz1 = model.dim_gnsf_nz1;
        ocp.dims.gnsf_nout = model.dim_gnsf_nout;
        ocp.dims.gnsf_ny = model.dim_gnsf_ny;
        ocp.dims.gnsf_nuhat = model.dim_gnsf_nuhat;
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

    ocp.cost.cost_ext_fun_type = model.cost_ext_fun_type;
    if strcmp(model.cost_ext_fun_type, 'generic')
        ocp.cost.cost_source_ext_cost = model.cost_source_ext_cost;
        ocp.cost.cost_function_ext_cost = model.cost_function_ext_cost;
    end
    ocp.cost.cost_ext_fun_type_0 = model.cost_ext_fun_type_0;
    if strcmp(model.cost_ext_fun_type_0, 'generic')
        ocp.cost.cost_source_ext_cost_0 = model.cost_source_ext_cost_0;
        ocp.cost.cost_function_ext_cost_0 = model.cost_function_ext_cost_0;
    end
    ocp.cost.cost_ext_fun_type_e = model.cost_ext_fun_type_e;
    if strcmp(model.cost_ext_fun_type_e, 'generic')
        ocp.cost.cost_source_ext_cost_e = model.cost_source_ext_cost_e;
        ocp.cost.cost_function_ext_cost_e = model.cost_function_ext_cost_e;
    end

    ocp.constraints.constr_type = upper(model.constr_type);
    ocp.constraints.constr_type_0 = upper(model.constr_type_0);
    ocp.constraints.constr_type_e = upper(model.constr_type_e);

    % parameters
    if model.dim_np > 0
        if isempty(obj.opts_struct.parameter_values)
            warning(['opts_struct.parameter_values are not set.', ...
                        10 'Using zeros(np,1) by default.' 10 'You can update them later using the solver object.']);
            ocp.parameter_values = zeros(model.dim_np,1);
        else
            ocp.parameter_values = obj.opts_struct.parameter_values(:);
        end
    end

    %% constraints
    % initial
    if isfield(model, 'constr_lbx_0')
        ocp.constraints.lbx_0 = model.constr_lbx_0;
    elseif ocp.dims.nbx_0 > 0
        warning('missing: constr_lbx_0, using zeros of appropriate dimension.');
        ocp.constraints.lbx_0 = zeros(ocp.dims.nbx_0, 1);
    end

    if isfield(model, 'constr_ubx_0')
        ocp.constraints.ubx_0 = model.constr_ubx_0;
    elseif ocp.dims.nbx_0 > 0
        warning('missing: constr_ubx_0, using zeros of appropriate dimension.');
        ocp.constraints.ubx_0 = zeros(ocp.dims.nbx_0, 1);
    end
    if isfield(model, 'constr_Jbx_0')
        ocp.constraints.idxbx_0 = J_to_idx( model.constr_Jbx_0 );
    end
    ocp.constraints.idxbxe_0 = model.constr_idxbxe_0;


    % path
    if ocp.dims.nbx > 0
        ocp.constraints.idxbx = J_to_idx( model.constr_Jbx );
        ocp.constraints.lbx = model.constr_lbx;
        ocp.constraints.ubx = model.constr_ubx;
    end

    if ocp.dims.nbu > 0
        ocp.constraints.idxbu = J_to_idx( model.constr_Jbu );
        ocp.constraints.lbu = model.constr_lbu;
        ocp.constraints.ubu = model.constr_ubu;
    end

    if ocp.dims.ng > 0
        ocp.constraints.C = model.constr_C;
        ocp.constraints.D = model.constr_D;
        ocp.constraints.lg = model.constr_lg;
        ocp.constraints.ug = model.constr_ug;
    end

    if ocp.dims.nh > 0
        ocp.constraints.lh = model.constr_lh;
        ocp.constraints.uh = model.constr_uh;
    end

    if ocp.dims.nsbx > 0
        ocp.constraints.idxsbx = J_to_idx_slack(model.constr_Jsbx);
        if isfield(model, 'constr_lsbx')
            ocp.constraints.lsbx = model.constr_lsbx;
        else
            ocp.constraints.lsbx = zeros(ocp.dims.nsbx, 1);
        end
        if isfield(model, 'constr_usbx')
            ocp.constraints.usbx = model.constr_usbx;
        else
            ocp.constraints.usbx = zeros(ocp.dims.nsbx, 1);
        end
    end


    if ocp.dims.nsbu > 0
        ocp.constraints.idxsbu = J_to_idx_slack(model.constr_Jsbu);
        if isfield(model, 'constr_lsbu')
            ocp.constraints.lsbu = model.constr_lsbu;
        else
            ocp.constraints.lsbu = zeros(ocp.dims.nsbu, 1);
        end
        if isfield(model, 'constr_usbu')
            ocp.constraints.usbu = model.constr_usbu;
        else
            ocp.constraints.usbu = zeros(ocp.dims.nsbu, 1);
        end
    end

    if ocp.dims.nsh > 0
        ocp.constraints.idxsh = J_to_idx_slack(model.constr_Jsh);
        if isfield(model, 'constr_lsh')
            ocp.constraints.lsh = model.constr_lsh;
        else
            ocp.constraints.lsh = zeros(ocp.dims.nsh, 1);
        end
        if isfield(model, 'constr_ush')
            ocp.constraints.ush = model.constr_ush;
        else
            ocp.constraints.ush = zeros(ocp.dims.nsh, 1);
        end
    end

    if ocp.dims.nsg > 0
        ocp.constraints.idxsg = J_to_idx_slack(model.constr_Jsg);
        if isfield(model, 'constr_lsg')
            ocp.constraints.lsg = model.constr_lsg;
        else
            ocp.constraints.lsg = zeros(ocp.dims.nsg, 1);
        end
        if isfield(model, 'constr_usg')
            ocp.constraints.usg = model.constr_usg;
        else
            ocp.constraints.usg = zeros(ocp.dims.nsg, 1);
        end
    end

    % terminal
    if ocp.dims.nbx_e > 0
        ocp.constraints.idxbx_e = J_to_idx( model.constr_Jbx_e );
        ocp.constraints.lbx_e = model.constr_lbx_e;
        ocp.constraints.ubx_e = model.constr_ubx_e;
    end

    if ocp.dims.ng_e > 0
        ocp.constraints.C_e = model.constr_C_e;
        ocp.constraints.lg_e = model.constr_lg_e;
        ocp.constraints.ug_e = model.constr_ug_e;
    end

    if ocp.dims.nh_e > 0
        ocp.constraints.lh_e = model.constr_lh_e;
        ocp.constraints.uh_e = model.constr_uh_e;
    end

    if ocp.dims.nh_0 > 0
        ocp.constraints.lh_0 = model.constr_lh_0;
        ocp.constraints.uh_0 = model.constr_uh_0;
    end

    if ocp.dims.nsbx_e > 0
        ocp.constraints.idxsbx_e = J_to_idx_slack(model.constr_Jsbx_e);
        if isfield(model, 'constr_lsbx_e')
            ocp.constraints.lsbx_e = model.constr_lsbx_e;
        else
            ocp.constraints.lsbx_e = zeros(ocp.dims.nsbx_e, 1);
        end
        if isfield(model, 'constr_usbx_e')
            ocp.constraints.usbx_e = model.constr_usbx_e;
        else
            ocp.constraints.usbx_e = zeros(ocp.dims.nsbx_e, 1);
        end
    end

    if ocp.dims.nsg_e > 0
        ocp.constraints.idxsg_e = J_to_idx_slack(model.constr_Jsg_e);
        if isfield(model, 'constr_lsg_e')
            ocp.constraints.lsg_e = model.constr_lsg_e;
        else
            ocp.constraints.lsg_e = zeros(ocp.dims.nsg_e, 1);
        end
        if isfield(model, 'constr_usg_e')
            ocp.constraints.usg_e = model.constr_usg_e;
        else
            ocp.constraints.usg_e = zeros(ocp.dims.nsg_e, 1);
        end
    end


    if ocp.dims.nsh_e > 0
        ocp.constraints.idxsh_e = J_to_idx_slack(model.constr_Jsh_e);
        if isfield(model, 'constr_lsh_e')
            ocp.constraints.lsh_e = model.constr_lsh_e;
        else
            ocp.constraints.lsh_e = zeros(ocp.dims.nsh_e, 1);
        end
        if isfield(model, 'constr_ush_e')
            ocp.constraints.ush_e = model.constr_ush_e;
        else
            ocp.constraints.ush_e = zeros(ocp.dims.nsh_e, 1);
        end
    end


    if ocp.dims.nsh_0 > 0
        ocp.constraints.idxsh_0 = J_to_idx_slack(model.constr_Jsh_0);
        if isfield(model, 'constr_lsh_0')
            ocp.constraints.lsh_0 = model.constr_lsh_0;
        else
            ocp.constraints.lsh_0 = zeros(ocp.dims.nsh_0, 1);
        end
        if isfield(model, 'constr_ush_0')
            ocp.constraints.ush_0 = model.constr_ush_0;
        else
            ocp.constraints.ush_0 = zeros(ocp.dims.nsh_0, 1);
        end
    end


    %% Cost
    if strcmp(model.cost_type, 'linear_ls')
        ocp.cost.Vu = model.cost_Vu;
        ocp.cost.Vx = model.cost_Vx;
        if isfield(model, 'cost_Vz')
            ocp.cost.Vz = model.cost_Vz;
        end
    end

    if strcmp(model.cost_type, 'nonlinear_ls') || strcmp(model.cost_type, 'linear_ls')
        ocp.cost.W = model.cost_W;
        if isfield(model, 'cost_y_ref')
            ocp.cost.yref = model.cost_y_ref;
        else
			warning(['cost_y_ref not defined for ocp json.' 10 'Using zeros(ny,1) by default.']);
            ocp.cost.yref = zeros(model.dim_ny,1);
        end
    end

    if strcmp(model.cost_type_0, 'linear_ls')
        ocp.cost.Vu_0 = model.cost_Vu_0;
        ocp.cost.Vx_0 = model.cost_Vx_0;
        if isfield(model, 'cost_Vz_0')
            ocp.cost.Vz_0 = model.cost_Vz_0;
        end
    end
    if strcmp(model.cost_type_0, 'nonlinear_ls') || strcmp(model.cost_type_0, 'linear_ls')
        ocp.cost.W_0 = model.cost_W_0;

        if isfield(model, 'cost_y_ref_0')
            ocp.cost.yref_0 = model.cost_y_ref_0;
        else
			warning(['cost_y_ref_0 not defined for ocp json.' 10 'Using zeros(ny_0,1) by default.']);
            ocp.cost.yref_0 = zeros(model.dim_ny_0,1);
        end
    end

    if isfield(model, 'cost_Vx_e')
        ocp.cost.Vx_e = model.cost_Vx_e;
    end

    if strcmp(model.cost_type_e, 'nonlinear_ls') || strcmp(model.cost_type_e, 'linear_ls')
        if isfield(model, 'cost_W_e')
            ocp.cost.W_e = model.cost_W_e;

        end
        if isfield(model, 'cost_y_ref_e')
            ocp.cost.yref_e = model.cost_y_ref_e;
        else
			warning(['cost_y_ref_e not defined for ocp json.' 10 'Using zeros(ny_e,1) by default.']);
            ocp.cost.yref_e = zeros(model.dim_ny_e,1);
        end
    end

    if isfield(model, 'cost_Zl')
        ocp.cost.Zl = diag(model.cost_Zl);
    end
    if isfield(model, 'cost_Zu')
        ocp.cost.Zu = diag(model.cost_Zu);
    end
    if isfield(model, 'cost_zl')
        ocp.cost.zl = model.cost_zl;
    end
    if isfield(model, 'cost_zu')
        ocp.cost.zu = model.cost_zu;
    end

    if isfield(model, 'cost_Zl_0')
        ocp.cost.Zl_0 = diag(model.cost_Zl_0);
    end
    if isfield(model, 'cost_Zu_0')
        ocp.cost.Zu_0 = diag(model.cost_Zu_0);
    end
    if isfield(model, 'cost_zl_0')
        ocp.cost.zl_0 = model.cost_zl_0;
    end
    if isfield(model, 'cost_zu_0')
        ocp.cost.zu_0 = model.cost_zu_0;
    end


    if isfield(model, 'cost_Zl_e')
        ocp.cost.Zl_e = diag(model.cost_Zl_e);
    end
    if isfield(model, 'cost_Zu_e')
        ocp.cost.Zu_e = diag(model.cost_Zu_e);
    end
    if isfield(model, 'cost_zl_e')
        ocp.cost.zl_e = model.cost_zl_e;
    end
    if isfield(model, 'cost_zu_e')
        ocp.cost.zu_e = model.cost_zu_e;
    end

    %% dynamics
    if strcmp(obj.opts_struct.sim_method, 'erk')
        ocp.model.f_expl_expr = model.dyn_expr_f;
    elseif strcmp(obj.opts_struct.sim_method, 'irk')
        if strcmp(model.dyn_ext_fun_type, 'casadi')
            ocp.model.f_impl_expr = model.dyn_expr_f;
        elseif strcmp(model.dyn_ext_fun_type, 'generic')
            ocp.model.dyn_generic_source = model.dyn_generic_source;
        end
    elseif strcmp(obj.opts_struct.sim_method, 'discrete')
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
    elseif strcmp(obj.opts_struct.sim_method, 'irk_gnsf')
        ocp.model.gnsf.A = model.dyn_gnsf_A;
        ocp.model.gnsf.B = model.dyn_gnsf_B;
        ocp.model.gnsf.C = model.dyn_gnsf_C;
        ocp.model.gnsf.E = model.dyn_gnsf_E;
        ocp.model.gnsf.c = model.dyn_gnsf_c;
        ocp.model.gnsf.A_LO = model.dyn_gnsf_A_LO;
        ocp.model.gnsf.B_LO = model.dyn_gnsf_B_LO;
        ocp.model.gnsf.E_LO = model.dyn_gnsf_E_LO;
        ocp.model.gnsf.c_LO = model.dyn_gnsf_c_LO;

        ocp.model.gnsf.L_x = model.dyn_gnsf_L_x;
        ocp.model.gnsf.L_u = model.dyn_gnsf_L_u;
        ocp.model.gnsf.L_xdot = model.dyn_gnsf_L_xdot;
        ocp.model.gnsf.L_z = model.dyn_gnsf_L_z;

        ocp.model.gnsf.expr_phi = model.dyn_gnsf_expr_phi;
        ocp.model.gnsf.expr_f_lo = model.dyn_gnsf_expr_f_lo;

        ocp.model.gnsf.ipiv_x = model.dyn_gnsf_ipiv_x;
        ocp.model.gnsf.idx_perm_x = model.dyn_gnsf_idx_perm_x;
        ocp.model.gnsf.ipiv_z = model.dyn_gnsf_ipiv_z;
        ocp.model.gnsf.idx_perm_z = model.dyn_gnsf_idx_perm_z;
        ocp.model.gnsf.ipiv_f = model.dyn_gnsf_ipiv_f;
        ocp.model.gnsf.idx_perm_f = model.dyn_gnsf_idx_perm_f;

        ocp.model.gnsf.nontrivial_f_LO = model.dyn_gnsf_nontrivial_f_LO;
        ocp.model.gnsf.purely_linear = model.dyn_gnsf_purely_linear;

        ocp.model.gnsf.y = model.sym_gnsf_y;
        ocp.model.gnsf.uhat = model.sym_gnsf_uhat;
    else
        error(['integrator ', obj.opts_struct.sim_method, ' not support for templating backend.'])
    end

    ocp.model.x = model.sym_x;
    ocp.model.u = model.sym_u;
    ocp.model.z = model.sym_z;
    ocp.model.xdot = model.sym_xdot;
    ocp.model.p = model.sym_p;

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
