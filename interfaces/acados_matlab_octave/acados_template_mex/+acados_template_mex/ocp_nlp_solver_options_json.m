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
        sim_method_jac_reuse
        nlp_solver_max_iter
        qp_solver_cond_N
        nlp_solver_tol_stat
        nlp_solver_tol_eq
        nlp_solver_tol_ineq
        nlp_solver_tol_comp
        nlp_solver_step_length
        rti_phase
        qp_solver_iter_max
        print_level
        initialize_t_slacks
        levenberg_marquardt
        regularize_method
        exact_hess_cost
        exact_hess_constr
        exact_hess_dyn
        ext_cost_num_hess

    end
    methods
        function obj = ocp_nlp_solver_options_json()
            obj.qp_solver       = 'PARTIAL_CONDENSING_HPIPM';
            obj.hessian_approx  = 'GAUSS_NEWTON';
            obj.integrator_type = 'ERK';
            obj.tf = [];
            obj.Tsim = [];
            obj.nlp_solver_type = 'SQP_RTI';
            obj.sim_method_num_steps = 1;
            obj.sim_method_num_stages = 2;
            obj.sim_method_newton_iter = 3;
            obj.sim_method_jac_reuse = 0;
            obj.nlp_solver_max_iter = 50;
            obj.qp_solver_cond_N = [];
            obj.nlp_solver_step_length = 1.0;
            obj.rti_phase = 0;
            obj.qp_solver_iter_max = [];
            obj.print_level = 0;
            obj.time_steps = [];
            obj.initialize_t_slacks = 0;
            obj.levenberg_marquardt = 0.0;
            obj.regularize_method = 'NO_REGULARIZE';
            obj.exact_hess_cost = 1;
            obj.exact_hess_constr = 1;
            obj.exact_hess_dyn = 1;
            obj.ext_cost_num_hess = 0;

        end
    end
end
