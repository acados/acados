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

classdef acados_ocp_opts < handle
	


	properties
		opts_struct
	end % properties



	methods
		

		function obj = acados_ocp_opts()
			% model stuct
			obj.opts_struct = struct;
			% default values
			obj.opts_struct.compile_mex = 'true';
			obj.opts_struct.codgen_model = 'true';
			obj.opts_struct.param_scheme = 'multiple_shooting_unif_grid';
			obj.opts_struct.param_scheme_N = 10;
			obj.opts_struct.nlp_solver = 'sqp';
			obj.opts_struct.nlp_solver_exact_hessian = 'false';
			obj.opts_struct.qp_solver = 'partial_condensing_hpipm';
			obj.opts_struct.sim_method = 'irk';
			obj.opts_struct.regularize_method = 'no_regularize';
			obj.opts_struct.gnsf_detect_struct = 'true';
			obj.opts_struct.output_dir = fullfile(pwd, 'build');
		end


		function obj = set(obj, field, value)
			if (strcmp(field, 'compile_mex'))
				obj.opts_struct.compile_mex = value;
			elseif (strcmp(field, 'codgen_model'))
				obj.opts_struct.codgen_model = value;
			elseif (strcmp(field, 'param_scheme'))
				obj.opts_struct.param_scheme = value;
			elseif (strcmp(field, 'param_scheme_N'))
				obj.opts_struct.param_scheme_N = value;
			elseif (strcmp(field, 'param_scheme_shooting_nodes'))
				obj.opts_struct.param_scheme_shooting_nodes = value;
			elseif (strcmp(field, 'nlp_solver'))
				obj.opts_struct.nlp_solver = value;
			elseif (strcmp(field, 'nlp_solver_exact_hessian'))
				obj.opts_struct.nlp_solver_exact_hessian = value;
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
			elseif (strcmp(field, 'qp_solver'))
				obj.opts_struct.qp_solver = value;
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
			elseif (strcmp(field, 'sim_method_num_stages'))
				obj.opts_struct.sim_method_num_stages = value;
			elseif (strcmp(field, 'sim_method_num_steps'))
				obj.opts_struct.sim_method_num_steps = value;
			elseif (strcmp(field, 'gnsf_detect_struct'))
				obj.opts_struct.gnsf_detect_struct = value;
			elseif (strcmp(field, 'regularize_method'))
				obj.opts_struct.regularize_method = value;
			elseif (strcmp(field, 'output_dir'))
				obj.opts_struct.output_dir = value;
			else
				disp(['acados_ocp_opts: set: wrong field: ', field]);
				keyboard;
			end
		end


	end % methods



end % class
