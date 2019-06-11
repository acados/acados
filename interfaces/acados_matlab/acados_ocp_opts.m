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
			elseif (strcmp(field, 'regularize_method'))
				obj.opts_struct.regularize_method = value;
			else
				disp(['acados_ocp_opts: set: wrong field: ', field]);
				keyboard;
			end
		end


	end % methods



end % class
