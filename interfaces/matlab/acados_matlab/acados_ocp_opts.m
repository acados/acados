classdef acados_ocp_opts < handle
	


	properties
		codgen_model
		nlp_solver
		qp_solver
		qp_solver_N_pcond
		sim_solver
		sim_solver_num_stages
		sim_solver_num_steps
		opts_struct
	end % properties



	methods
		

		function obj = acados_ocp_opts()
			obj.codgen_model = 'true';
			nlp_solver = 0;
			qp_solver = 0;
			qp_solver_N_pcond = 0;
			obj.sim_solver = 0;
			obj.sim_solver_num_stages = 0;
			obj.sim_solver_num_steps = 0;
			obj.opts_struct = struct;
			obj.opts_struct.codgen_model = obj.codgen_model;
		end


		function obj = set(obj, field, value)
			if (strcmp(field, 'codgen_model'))
				obj.codgen_model = value;
				obj.opts_struct.codgen_model = value;
			elseif (strcmp(field, 'nlp_solver'))
				obj.nlp_solver = value;
				obj.opts_struct.nlp_solver = value;
			elseif (strcmp(field, 'qp_solver'))
				obj.qp_solver = value;
				obj.opts_struct.qp_solver = value;
			elseif (strcmp(field, 'qp_solver_N_pcond'))
				obj.qp_solver_N_pcond = value;
				obj.opts_struct.qp_solver_N_pcond = value;
			elseif (strcmp(field, 'sim_solver'))
				obj.sim_solver = value;
				obj.opts_struct.sim_solver = value;
			elseif (strcmp(field, 'sim_solver_num_stages'))
				obj.sim_solver_num_stages = value;
				obj.opts_struct.sim_solver_num_stages = value;
			elseif (strcmp(field, 'sim_solver_num_steps'))
				obj.sim_solver_num_steps = value;
				obj.opts_struct.sim_solver_num_steps = value;
			else
				disp('acados_ocp_opts: set: wrong field');
			end
		end


	end % methods



end % class
