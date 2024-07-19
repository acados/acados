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
        nlp_solver_type        %  NLP solver
        nlp_solver_step_length
        nlp_solver_tol_stat
        nlp_solver_tol_eq
        nlp_solver_tol_ineq
        nlp_solver_tol_comp
        nlp_solver_max_iter
        nlp_solver_ext_qp_res
        nlp_solver_warm_start_first_qp
        nlp_solver_tol_min_step_norm
        globalization
        levenberg_marquardt
        collocation_type
        sim_method_num_steps   %  number of steps in integrator
        sim_method_num_stages  %  size of butcher tableau
        sim_method_newton_iter
        sim_method_newton_tol
        sim_method_jac_reuse
        time_steps
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
        rti_log_residuals
        print_level
        cost_discretization
        regularize_method
        reg_epsilon
        exact_hess_cost
        exact_hess_dyn
        exact_hess_constr
        fixed_hess
        ext_cost_num_hess
        alpha_min
        alpha_reduction
        line_search_use_sufficient_descent
        globalization_use_SOC
        full_step_dual
        eps_sufficient_descent
        hpipm_mode
        with_solution_sens_wrt_params
        with_value_sens_wrt_params
        as_rti_iter
        as_rti_level
        with_adaptive_levenberg_marquardt
        adaptive_levenberg_marquardt_lam
        adaptive_levenberg_marquardt_mu_min
        adaptive_levenberg_marquardt_mu0
        log_primal_step_norm
        eval_residual_at_max_iter

        ext_fun_compile_flags
        model_external_shared_lib_dir
        model_external_shared_lib_name
        custom_update_filename
        custom_update_header_filename
        custom_templates
        custom_update_copy
        num_threads_in_batch_solve

    end
    methods
        function obj = AcadosOcpOptions()
            obj.hessian_approx = 'GAUSS_NEWTON';
            obj.integrator_type = 'ERK';
            obj.tf = [];
            obj.nlp_solver_type = 'SQP_RTI';
            obj.nlp_solver_step_length = 1.0;
            obj.nlp_solver_tol_stat = 1e-6;
            obj.nlp_solver_tol_eq = 1e-6;
            obj.nlp_solver_tol_ineq = 1e-6;
            obj.nlp_solver_tol_comp = 1e-6;
            obj.nlp_solver_tol_min_step_norm = 1e-12;
            obj.nlp_solver_max_iter = 50;
            obj.nlp_solver_ext_qp_res = 0;
            obj.nlp_solver_warm_start_first_qp = false;
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
            obj.qp_solver_cond_ric_alg = 0;
            obj.qp_solver_ric_alg = 0;
            obj.rti_log_residuals = 0;
            obj.print_level = 0;
            obj.cost_discretization = 'EULER';
            obj.regularize_method = 'NO_REGULARIZE';
            obj.reg_epsilon = 1e-4;
            obj.exact_hess_cost = 1;
            obj.exact_hess_dyn = 1;
            obj.exact_hess_constr = 1;
            obj.fixed_hess = 0;
            obj.ext_cost_num_hess = 0;
            obj.alpha_min = 0.05;
            obj.alpha_reduction = 0.7;
            obj.line_search_use_sufficient_descent = 0;
            obj.globalization_use_SOC = 0;
            obj.full_step_dual = 0;
            obj.eps_sufficient_descent = 1e-4;
            obj.hpipm_mode = 'BALANCE';
            obj.with_solution_sens_wrt_params = 0;
            obj.with_value_sens_wrt_params = 0;
            obj.as_rti_iter = 1;
            obj.as_rti_level = 4;
            obj.with_adaptive_levenberg_marquardt = 0;
            obj.adaptive_levenberg_marquardt_lam = 5.0;
            obj.adaptive_levenberg_marquardt_mu_min = 1e-16;
            obj.adaptive_levenberg_marquardt_mu0 = 1e-3;
            obj.log_primal_step_norm = 0;
            obj.eval_residual_at_max_iter = 0;

            obj.ext_fun_compile_flags = '-O2';
            obj.model_external_shared_lib_dir = [];
            obj.model_external_shared_lib_name = [];
            obj.custom_update_filename = '';
            obj.custom_update_header_filename = '';
            obj.custom_templates = [];
            obj.custom_update_copy = true;
            obj.num_threads_in_batch_solve = 1;
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
    end
end
