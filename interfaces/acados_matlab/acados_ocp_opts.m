classdef acados_ocp_opts < handle
	


	properties
		compile_mex
		codgen_model
		param_scheme
		param_scheme_N
		nlp_solver
		qp_solver
		qp_solver_N_pcond
		sim_method
		sim_method_num_stages
		sim_method_num_steps
		opts_struct
	end % properties



	methods
		

		function obj = acados_ocp_opts()
			% default values
			obj.compile_mex = 'true';
			obj.codgen_model = 'true';
			obj.param_scheme = 'multiple_shooting_unif_grid';
			obj.param_scheme_N = 10;
			obj.nlp_solver = 'sqp';
			obj.qp_solver = 'qp_solver';
			obj.sim_method = 'irk';
			% model stuct
			obj.opts_struct = struct;
			% initialize model stuct
			obj.opts_struct.compile_mex = obj.compile_mex;
			obj.opts_struct.codgen_model = obj.codgen_model;
			obj.opts_struct.param_scheme = obj.param_scheme;
			obj.opts_struct.param_scheme_N = obj.param_scheme_N;
			obj.opts_struct.qp_solver = obj.qp_solver;
			obj.opts_struct.nlp_solver = obj.nlp_solver;
			obj.opts_struct.sim_method = obj.sim_method;
		end


		function obj = set(obj, field, value)
			if (strcmp(field, 'compile_mex'))
				obj.compile_mex = value;
				obj.opts_struct.compile_mex = value;
			elseif (strcmp(field, 'codgen_model'))
				obj.codgen_model = value;
				obj.opts_struct.codgen_model = value;
			elseif (strcmp(field, 'param_scheme'))
				obj.param_scheme = value;
				obj.opts_struct.param_scheme = value;
			elseif (strcmp(field, 'param_scheme_N'))
				obj.param_scheme_N = value;
				obj.opts_struct.param_scheme_N = value;
			elseif (strcmp(field, 'nlp_solver'))
				obj.nlp_solver = value;
				obj.opts_struct.nlp_solver = value;
			elseif (strcmp(field, 'qp_solver'))
				obj.qp_solver = value;
				obj.opts_struct.qp_solver = value;
			elseif (strcmp(field, 'qp_solver_N_pcond'))
				obj.qp_solver_N_pcond = value;
				obj.opts_struct.qp_solver_N_pcond = value;
			elseif (strcmp(field, 'sim_method'))
				obj.sim_method = value;
				obj.opts_struct.sim_method = value;
			elseif (strcmp(field, 'sim_method_num_stages'))
				obj.sim_method_num_stages = value;
				obj.opts_struct.sim_method_num_stages = value;
			elseif (strcmp(field, 'sim_method_num_steps'))
				obj.sim_method_num_steps = value;
				obj.opts_struct.sim_method_num_steps = value;
			else
				disp('acados_ocp_opts: set: wrong field');
			end
		end


	end % methods



end % class
