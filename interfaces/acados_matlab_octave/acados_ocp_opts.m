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

classdef acados_ocp_opts < handle


    properties
        opts_struct
    end % properties



    methods

        function obj = acados_ocp_opts()
            % model stuct
            obj.opts_struct = struct;
            % default values
            obj.opts_struct.compile_interface = 'auto'; % auto, true, false
            obj.opts_struct.codgen_model = 'true';
            obj.opts_struct.compile_model = 'true';
            obj.opts_struct.param_scheme_N = 10;
            % set one of the following for nonuniform grid
            obj.opts_struct.shooting_nodes = [];
            obj.opts_struct.time_steps = [];
            obj.opts_struct.parameter_values = [];

            obj.opts_struct.nlp_solver = 'sqp';
            obj.opts_struct.nlp_solver_exact_hessian = 'false';
            obj.opts_struct.nlp_solver_max_iter = 100;
            obj.opts_struct.nlp_solver_tol_stat = 1e-6;
            obj.opts_struct.nlp_solver_tol_eq = 1e-6;
            obj.opts_struct.nlp_solver_tol_ineq = 1e-6;
            obj.opts_struct.nlp_solver_tol_comp = 1e-6;
            obj.opts_struct.nlp_solver_ext_qp_res = 0; % compute QP residuals at each NLP iteration
            obj.opts_struct.nlp_solver_step_length = 1.0; % fixed step length in SQP algorithm
            obj.opts_struct.rti_phase = 0; % RTI phase: (1) preparation, (2) feedback, (0) both
            obj.opts_struct.qp_solver = 'partial_condensing_hpipm';
            % globalization
            obj.opts_struct.globalization = 'fixed_step';
            obj.opts_struct.alpha_min = 0.05;
            obj.opts_struct.alpha_reduction = 0.7;
            obj.opts_struct.line_search_use_sufficient_descent = 0;
            obj.opts_struct.globalization_use_SOC = 0;
            obj.opts_struct.full_step_dual = 0;
            obj.opts_struct.eps_sufficient_descent = 1e-4;

            obj.opts_struct.qp_solver_iter_max = 50;
            % obj.opts_struct.qp_solver_cond_N = 5; % New horizon after partial condensing
            obj.opts_struct.qp_solver_cond_ric_alg = 1; % 0: dont factorize hessian in the condensing; 1: factorize
            obj.opts_struct.qp_solver_ric_alg = 1; % HPIPM specific
            obj.opts_struct.qp_solver_warm_start = 0;
                    % 0 no warm start; 1 warm start primal variables; 2 warm start primal and dual variables
            obj.opts_struct.warm_start_first_qp = 0;
                    % 0 no warm start in first sqp iter - 1 warm start even in first sqp iter
            obj.opts_struct.sim_method = 'irk'; % erk; irk; irk_gnsf
            obj.opts_struct.collocation_type = 'gauss_legendre';
            obj.opts_struct.sim_method_num_stages = 4;
            obj.opts_struct.sim_method_num_steps = 1;
            obj.opts_struct.sim_method_newton_iter = 3;
            obj.opts_struct.sim_method_jac_reuse = 0;
            obj.opts_struct.gnsf_detect_struct = 'true';
            obj.opts_struct.regularize_method = 'no_regularize';
            obj.opts_struct.reg_epsilon = 1e-4;
            obj.opts_struct.print_level = 0;
            obj.opts_struct.levenberg_marquardt = 0.0;
            % 0 or 1, only used if nlp_solver_exact_hessian
            obj.opts_struct.exact_hess_dyn = 1;
            obj.opts_struct.exact_hess_cost = 1;
            obj.opts_struct.exact_hess_constr = 1;
            obj.opts_struct.fixed_hess = 0;
            obj.opts_struct.ext_fun_compile_flags = '-O2';

            obj.opts_struct.output_dir = fullfile(pwd, 'build');
            % if ismac()
            %     obj.opts_struct.output_dir = '/usr/local/lib';
            % end
        end


        function obj = set(obj, field, value)
            % convert Matlab strings to char arrays
            if isstring(value)
                value = char(value);
            end

            if (strcmp(field, 'compile_interface'))
                obj.opts_struct.compile_interface = value;
            elseif (strcmp(field, 'codgen_model'))
                obj.opts_struct.codgen_model = value;
            elseif (strcmp(field, 'compile_model'))
                obj.opts_struct.compile_model = value;
            elseif (strcmp(field, 'param_scheme'))
                warning(['param_scheme: option is outdated! Uniform discretization with T/N is default!\n',...
                         'Set opts.shooting_nodes or opts.time_steps for nonuniform discretizations.'])
            elseif (strcmp(field, 'param_scheme_N'))
                obj.opts_struct.param_scheme_N = value;
            elseif (any(strcmp(field, {'param_scheme_shooting_nodes','shooting_nodes'})))
                obj.opts_struct.shooting_nodes = value;
            elseif (strcmp(field, 'time_steps'))
                obj.opts_struct.time_steps = value;
            elseif (strcmp(field, 'nlp_solver'))
                obj.opts_struct.nlp_solver = value;
            elseif (strcmp(field, 'nlp_solver_exact_hessian'))
                obj.opts_struct.nlp_solver_exact_hessian = value;
            % hessian approx
            elseif (strcmp(field, 'exact_hess_dyn'))
                obj.opts_struct.exact_hess_dyn = value;
            elseif (strcmp(field, 'exact_hess_cost'))
                obj.opts_struct.exact_hess_cost = value;
            elseif (strcmp(field, 'exact_hess_constr'))
                obj.opts_struct.exact_hess_constr = value;
            elseif (strcmp(field, 'fixed_hess'))
                obj.opts_struct.fixed_hess = value;
            elseif (strcmp(field, 'nlp_solver_max_iter'))
                obj.opts_struct.nlp_solver_max_iter = value;
            elseif (strcmp(field, 'nlp_solver_tol_stat'))
                obj.opts_struct.nlp_solver_tol_stat = value;
            elseif (strcmp(field, 'nlp_solver_tol_eq'))
                obj.opts_struct.nlp_solver_tol_eq = value;
            elseif (strcmp(field, 'nlp_solver_tol_ineq'))
                obj.opts_struct.nlp_solver_tol_ineq = value;
            elseif (strcmp(field, 'nlp_solver_tol_comp'))
                obj.opts_struct.nlp_solver_tol_comp = value;
            elseif (strcmp(field, 'nlp_solver_ext_qp_res'))
                obj.opts_struct.nlp_solver_ext_qp_res = value;
            elseif (strcmp(field, 'nlp_solver_step_length'))
                obj.opts_struct.nlp_solver_step_length = value;
            elseif (strcmp(field, 'rti_phase'))
                obj.opts_struct.rti_phase = value;
            elseif (strcmp(field, 'nlp_solver_warm_start_first_qp'))
                obj.opts_struct.nlp_solver_warm_start_first_qp = value;
            elseif (strcmp(field, 'qp_solver'))
                obj.opts_struct.qp_solver = value;
            elseif (strcmp(field, 'qp_solver_tol_stat'))
				obj.opts_struct.qp_solver_tol_stat = value;
			elseif (strcmp(field, 'qp_solver_tol_eq'))
				obj.opts_struct.qp_solver_tol_eq = value;
			elseif (strcmp(field, 'qp_solver_tol_ineq'))
				obj.opts_struct.qp_solver_tol_ineq = value;
			elseif (strcmp(field, 'qp_solver_tol_comp'))
				obj.opts_struct.qp_solver_tol_comp = value;
            elseif (strcmp(field, 'qp_solver_iter_max'))
                obj.opts_struct.qp_solver_iter_max = value;
            elseif (strcmp(field, 'qp_solver_cond_N'))
                obj.opts_struct.qp_solver_cond_N = value;
            elseif (strcmp(field, 'qp_solver_cond_ric_alg'))
                obj.opts_struct.qp_solver_cond_ric_alg = value;
            elseif (strcmp(field, 'qp_solver_ric_alg'))
                obj.opts_struct.qp_solver_ric_alg = value;
            elseif (strcmp(field, 'qp_solver_warm_start'))
                obj.opts_struct.qp_solver_warm_start = value;
            elseif (strcmp(field, 'sim_method'))
                obj.opts_struct.sim_method = value;
            elseif (strcmp(field, 'collocation_type'))
                obj.opts_struct.collocation_type = value;
            elseif (strcmp(field, 'sim_method_num_stages'))
                obj.opts_struct.sim_method_num_stages = value;
            elseif (strcmp(field, 'sim_method_num_steps'))
                obj.opts_struct.sim_method_num_steps = value;
            elseif (strcmp(field, 'sim_method_newton_iter'))
                obj.opts_struct.sim_method_newton_iter = value;
            elseif (strcmp(field, 'sim_method_exact_z_output'))
                obj.opts_struct.sim_method_exact_z_output = value;
            elseif (strcmp(field, 'sim_method_jac_reuse'))
                obj.opts_struct.sim_method_jac_reuse = value;
            elseif (strcmp(field, 'gnsf_detect_struct'))
                obj.opts_struct.gnsf_detect_struct = value;
            elseif (strcmp(field, 'regularize_method'))
                obj.opts_struct.regularize_method = value;
            elseif (strcmp(field, 'reg_epsilon'))
                obj.opts_struct.reg_epsilon = value;
            elseif (strcmp(field, 'output_dir'))
                obj.opts_struct.output_dir = value;
            elseif (strcmp(field, 'print_level'))
                obj.opts_struct.print_level = value;
            elseif (strcmp(field, 'levenberg_marquardt'))
                obj.opts_struct.levenberg_marquardt = value;
            elseif (strcmp(field, 'alpha_min'))
                obj.opts_struct.alpha_min = value;
            elseif (strcmp(field, 'alpha_reduction'))
                obj.opts_struct.alpha_reduction = value;
            elseif (strcmp(field, 'line_search_use_sufficient_descent'))
                obj.opts_struct.line_search_use_sufficient_descent = value;
            elseif (strcmp(field, 'globalization_use_SOC'))
                obj.opts_struct.globalization_use_SOC = value;
            elseif (strcmp(field, 'full_step_dual'))
                obj.opts_struct.full_step_dual = value;
            elseif (strcmp(field, 'eps_sufficient_descent'))
                obj.opts_struct.eps_sufficient_descent = value;
            elseif (strcmp(field, 'globalization'))
                obj.opts_struct.globalization = value;
            elseif (strcmp(field, 'parameter_values'))
                obj.opts_struct.parameter_values = value;
            elseif (strcmp(field, 'ext_fun_compile_flags'))
                obj.opts_struct.ext_fun_compile_flags = value;
            elseif (strcmp(field, 'compile_mex'))
                disp(['Option compile_mex is not supported anymore,'...
                    'please use compile_interface instead or dont set the option.', ...
                    'options are: true, false, auto.']);
                keyboard
            else
                disp(['acados_ocp_opts: set: wrong field: ', field]);
                keyboard;
            end
        end

    end % methods


end % class
