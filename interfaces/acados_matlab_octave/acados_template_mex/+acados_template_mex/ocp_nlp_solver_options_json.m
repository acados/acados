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

classdef ocp_nlp_solver_options_json < handle
    properties
        qp_solver              %  qp solver to be used in the NLP solver
        hessian_approx         %  hessian approximation
        integrator_type        %  integrator type
        tf                     %  prediction horizon
        Tsim
        time_steps
        nlp_solver_type        %  NLP solver
        sim_method_num_steps   %  number of steps in integrator
        sim_method_num_stages  %  size of butcher tableau
        sim_method_newton_iter
        sim_method_newton_tol
        sim_method_jac_reuse
        nlp_solver_max_iter
        nlp_solver_tol_stat
        nlp_solver_tol_eq
        nlp_solver_tol_ineq
        nlp_solver_tol_comp
        nlp_solver_step_length
        nlp_solver_warm_start_first_qp
        rti_phase
        qp_solver_cond_N
        qp_solver_iter_max
        qp_solver_tol_stat
        reg_epsilon
        qp_solver_tol_eq
        qp_solver_tol_ineq
        qp_solver_tol_comp
        qp_solver_warm_start
        print_level
        initialize_t_slacks
        levenberg_marquardt
        regularize_method
        exact_hess_cost
        exact_hess_constr
        exact_hess_dyn
        fixed_hess
        ext_cost_num_hess
        alpha_min
        alpha_reduction
        globalization
        collocation_type
        line_search_use_sufficient_descent
        globalization_use_SOC
        full_step_dual
        eps_sufficient_descent
        hpipm_mode
        qp_solver_ric_alg
        qp_solver_cond_ric_alg
        nlp_solver_ext_qp_res
        ext_fun_compile_flags
        cost_discretization

    end
    methods
        function obj = ocp_nlp_solver_options_json()
            obj.qp_solver       = 'PARTIAL_CONDENSING_HPIPM';
            obj.hessian_approx  = 'GAUSS_NEWTON';
            obj.integrator_type = 'ERK';
            obj.collocation_type = 'GAUSS_LEGENDRE';
            obj.tf = [];
            obj.Tsim = [];
            obj.nlp_solver_type = 'SQP_RTI';
            obj.sim_method_num_steps = 1;
            obj.sim_method_num_stages = 2;
            obj.sim_method_newton_iter = 3;
            obj.sim_method_newton_tol = 0.0;
            obj.sim_method_jac_reuse = 0;
            obj.nlp_solver_max_iter = 50;
            obj.qp_solver_cond_N = [];
            obj.nlp_solver_step_length = 1.0;
            obj.rti_phase = 0;
            obj.qp_solver_iter_max = [];
            obj.reg_epsilon = 1e-4;
            obj.print_level = 0;
            obj.time_steps = [];
            obj.initialize_t_slacks = 0;
            obj.levenberg_marquardt = 0.0;
            obj.regularize_method = 'NO_REGULARIZE';
            obj.exact_hess_cost = 1;
            obj.exact_hess_constr = 1;
            obj.exact_hess_dyn = 1;
            obj.fixed_hess = 0;
            obj.ext_cost_num_hess = 0;
            obj.alpha_min = 0.05;
            obj.alpha_reduction = 0.7;
            obj.globalization = 'FIXED_STEP';
            obj.line_search_use_sufficient_descent = 0;
            obj.globalization_use_SOC = 0;
            obj.full_step_dual = 0;
            obj.qp_solver_cond_ric_alg = 0;
            obj.qp_solver_ric_alg = 0;
            obj.eps_sufficient_descent = 1e-4;
            obj.hpipm_mode = 'BALANCE';
            obj.nlp_solver_ext_qp_res = 0;
            obj.ext_fun_compile_flags = '-O2';
            obj.cost_discretization = 'EULER';

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
