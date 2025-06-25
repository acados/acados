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



classdef AcadosOcpOptions < handle
    properties
        hessian_approx         %  hessian approximation
        integrator_type        %  integrator type
        tf                     %  prediction horizon
        N_horizon

        nlp_solver_type        %  NLP solver
        nlp_solver_step_length % TODO: is deprecated, remove in future release
        nlp_solver_tol_stat
        nlp_solver_tol_eq
        nlp_solver_tol_ineq
        nlp_solver_tol_comp
        nlp_solver_max_iter
        nlp_solver_ext_qp_res
        nlp_solver_warm_start_first_qp
        nlp_solver_warm_start_first_qp_from_nlp
        nlp_solver_tol_min_step_norm
        globalization
        levenberg_marquardt
        collocation_type
        sim_method_num_steps   %  number of steps in integrator
        sim_method_num_stages  %  size of butcher tableau
        sim_method_newton_iter
        sim_method_newton_tol
        sim_method_jac_reuse
        sim_method_detect_gnsf
        time_steps
        shooting_nodes
        cost_scaling
        Tsim
        qp_solver              %  qp solver to be used in the NLP solver
        qp_solver_tol_stat
        qp_solver_tol_eq
        qp_solver_tol_ineq
        qp_solver_tol_comp
        qp_solver_iter_max
        qp_solver_cond_N
        qp_solver_cond_block_size
        qp_solver_warm_start
        qp_solver_cond_ric_alg
        qp_solver_ric_alg
        qp_solver_mu0
        qp_solver_t0_init
        tau_min
        rti_log_residuals
        rti_log_only_available_residuals
        print_level
        cost_discretization
        regularize_method
        reg_epsilon
        reg_max_cond_block
        reg_min_epsilon
        reg_adaptive_eps
        qpscaling_ub_max_abs_eig
        qpscaling_lb_norm_inf_grad_obj
        qpscaling_scale_objective
        qpscaling_scale_constraints
        exact_hess_cost
        exact_hess_dyn
        exact_hess_constr
        fixed_hess
        ext_cost_num_hess
        globalization_fixed_step_length
        globalization_alpha_min
        globalization_alpha_reduction
        globalization_line_search_use_sufficient_descent
        globalization_use_SOC
        globalization_full_step_dual
        globalization_eps_sufficient_descent
        globalization_funnel_init_increase_factor
        globalization_funnel_init_upper_bound
        globalization_funnel_sufficient_decrease_factor
        globalization_funnel_kappa
        globalization_funnel_fraction_switching_condition
        globalization_funnel_initial_penalty_parameter
        globalization_funnel_use_merit_fun_only


        search_direction_mode
        use_constraint_hessian_in_feas_qp
        allow_direction_mode_switch_to_nominal
        hpipm_mode
        with_solution_sens_wrt_params
        with_value_sens_wrt_params
        solution_sens_qp_t_lam_min
        as_rti_iter
        as_rti_level
        with_adaptive_levenberg_marquardt
        adaptive_levenberg_marquardt_lam
        adaptive_levenberg_marquardt_mu_min
        adaptive_levenberg_marquardt_mu0
        adaptive_levenberg_marquardt_obj_scalar
        log_primal_step_norm
        log_dual_step_norm
        store_iterates
        eval_residual_at_max_iter
        with_anderson_acceleration

        timeout_max_time
        timeout_heuristic

        ext_fun_compile_flags
        ext_fun_expand_dyn
        ext_fun_expand_cost
        ext_fun_expand_constr
        ext_fun_expand_precompute

        model_external_shared_lib_dir
        model_external_shared_lib_name
        custom_update_filename
        custom_update_header_filename
        custom_templates
        custom_update_copy
        with_batch_functionality

        compile_interface

    end
    methods
        function obj = AcadosOcpOptions()
            obj.hessian_approx = 'GAUSS_NEWTON';
            obj.integrator_type = 'ERK';
            obj.tf = [];
            obj.N_horizon = [];
            obj.nlp_solver_type = 'SQP';
            obj.globalization_fixed_step_length = 1.0;
            obj.nlp_solver_step_length = [];
            obj.nlp_solver_tol_stat = 1e-6;
            obj.nlp_solver_tol_eq = 1e-6;
            obj.nlp_solver_tol_ineq = 1e-6;
            obj.nlp_solver_tol_comp = 1e-6;
            obj.nlp_solver_tol_min_step_norm = [];
            obj.nlp_solver_max_iter = 100;
            obj.nlp_solver_ext_qp_res = 0;
            obj.nlp_solver_warm_start_first_qp = false;
            obj.nlp_solver_warm_start_first_qp_from_nlp = false;
            obj.globalization = 'FIXED_STEP';
            obj.levenberg_marquardt = 0.0;
            obj.collocation_type = 'GAUSS_LEGENDRE';
            obj.sim_method_num_stages = 4;
            obj.sim_method_num_steps = 1;
            obj.sim_method_newton_iter = 3;
            obj.sim_method_newton_tol = 0.0;
            obj.sim_method_jac_reuse = 0;
            obj.time_steps = [];
            obj.Tsim = [];
            obj.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
            obj.qp_solver_tol_stat = [];
            obj.qp_solver_tol_eq = [];
            obj.qp_solver_tol_ineq = [];
            obj.qp_solver_tol_comp = [];
            obj.qp_solver_iter_max = 50;
            obj.qp_solver_cond_N = [];
            obj.qp_solver_cond_block_size = [];
            obj.qp_solver_cond_ric_alg = 1;
            obj.qp_solver_ric_alg = 1;
            obj.qp_solver_mu0 = 0;
            obj.qp_solver_t0_init = 2;
            obj.tau_min = 0;
            obj.rti_log_residuals = 0;
            obj.rti_log_only_available_residuals = 0;
            obj.print_level = 0;
            obj.cost_discretization = 'EULER';
            obj.regularize_method = 'NO_REGULARIZE';
            obj.qpscaling_ub_max_abs_eig = 1e5;
            obj.qpscaling_lb_norm_inf_grad_obj = 1e-4;
            obj.qpscaling_scale_objective = 'NO_OBJECTIVE_SCALING';
            obj.qpscaling_scale_constraints = 'NO_CONSTRAINT_SCALING';
            obj.reg_epsilon = 1e-4;
            obj.reg_adaptive_eps = false;
            obj.reg_max_cond_block = 1e7;
            obj.reg_min_epsilon = 1e-8;
            obj.shooting_nodes = [];
            obj.cost_scaling = [];
            obj.exact_hess_cost = 1;
            obj.exact_hess_dyn = 1;
            obj.exact_hess_constr = 1;
            obj.fixed_hess = 0;
            obj.ext_cost_num_hess = 0;
            obj.globalization_alpha_min = [];
            obj.globalization_alpha_reduction = [];
            obj.globalization_line_search_use_sufficient_descent = 0;
            obj.globalization_use_SOC = 0;
            obj.globalization_full_step_dual = [];
            obj.globalization_eps_sufficient_descent = [];


            % funnel options
            obj.globalization_funnel_init_increase_factor = 15;
            obj.globalization_funnel_init_upper_bound = 1.0;
            obj.globalization_funnel_sufficient_decrease_factor = 0.9;
            obj.globalization_funnel_kappa = 0.9;
            obj.globalization_funnel_fraction_switching_condition = 1e-3;
            obj.globalization_funnel_initial_penalty_parameter = 1.0;
            obj.globalization_funnel_use_merit_fun_only = false;

            % SQP_WITH_FEASIBLE_QP options
            obj.search_direction_mode = 'NOMINAL_QP';
            obj.use_constraint_hessian_in_feas_qp = false;
            obj.allow_direction_mode_switch_to_nominal = true;

            obj.hpipm_mode = 'BALANCE';
            obj.with_solution_sens_wrt_params = 0;
            obj.with_value_sens_wrt_params = 0;
            obj.solution_sens_qp_t_lam_min = 1e-9;
            obj.as_rti_iter = 1;
            obj.as_rti_level = 4;
            obj.with_adaptive_levenberg_marquardt = 0;
            obj.adaptive_levenberg_marquardt_lam = 5.0;
            obj.adaptive_levenberg_marquardt_mu_min = 1e-16;
            obj.adaptive_levenberg_marquardt_mu0 = 1e-3;
            obj.adaptive_levenberg_marquardt_obj_scalar = 2.0;
            obj.log_primal_step_norm = 0;
            obj.log_dual_step_norm = 0;
            obj.store_iterates = false;
            obj.eval_residual_at_max_iter = [];
            obj.with_anderson_acceleration = 0;
            obj.timeout_max_time = 0.;
            obj.timeout_heuristic = 'ZERO';

            % check whether flags are provided by environment variable
            env_var = getenv("ACADOS_EXT_FUN_COMPILE_FLAGS");
            if isempty(env_var)
                obj.ext_fun_compile_flags = '-O2';
            else
                obj.ext_fun_compile_flags = env_var;
            end
            obj.ext_fun_expand_dyn = false;
            obj.ext_fun_expand_cost = false;
            obj.ext_fun_expand_constr = false;
            obj.ext_fun_expand_precompute = false;

            obj.model_external_shared_lib_dir = [];
            obj.model_external_shared_lib_name = [];
            obj.custom_update_filename = '';
            obj.custom_update_header_filename = '';
            obj.custom_templates = [];
            obj.custom_update_copy = true;
            obj.with_batch_functionality = false;

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

        function s = convert_to_struct_for_json_dump(self, N)
            s = self.struct();
            s = prepare_struct_for_json_dump(s, {'time_steps', 'shooting_nodes', 'cost_scaling', 'sim_method_num_stages', 'sim_method_num_steps', 'sim_method_jac_reuse', 'custom_templates'}, {});
        end
    end
end
