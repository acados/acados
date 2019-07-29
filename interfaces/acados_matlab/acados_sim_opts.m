classdef acados_sim_opts < handle
	


	properties
		opts_struct
	end % properties



	methods
		

		function obj = acados_sim_opts()
			obj.opts_struct = struct;
			obj.opts_struct.compile_mex = 'true';
			obj.opts_struct.codgen_model = 'true';
			obj.opts_struct.method = 'irk';
			obj.opts_struct.num_stages = 4;
			obj.opts_struct.num_steps = 1;
			obj.opts_struct.sens_forw = 'false';
			obj.opts_struct.sens_adj = 'false';
			obj.opts_struct.sens_hess = 'false';
			obj.opts_struct.gnsf_detect_struct = 'true';
			obj.opts_struct.output_dir = fullfile(pwd, 'build');
		end


		function obj = set(obj, field, value)
			if (strcmp(field, 'compile_mex'))
				obj.opts_struct.compile_mex = value;
			elseif (strcmp(field, 'codgen_model'))
				obj.opts_struct.codgen_model = value;
			elseif (strcmp(field, 'num_stages'))
				obj.opts_struct.num_stages = value;
			elseif (strcmp(field, 'num_steps'))
				obj.opts_struct.num_steps = value;
			elseif (strcmp(field, 'method'))
				obj.opts_struct.method = value;
			elseif (strcmp(field, 'sens_forw'))
				obj.opts_struct.sens_forw = value;
			elseif (strcmp(field, 'sens_adj'))
				obj.opts_struct.sens_adj = value;
			elseif (strcmp(field, 'sens_hess'))
				obj.opts_struct.sens_hess = value;
			elseif (strcmp(field, 'gnsf_detect_struct'))
				obj.opts_struct.gnsf_detect_struct = value;
			elseif (strcmp(field, 'output_dir'))
				obj.opts_struct.output_dir = value;
			else
				disp(['acados_sim_opts: set: wrong field: ', field]);
			end
		end


	end % methods



end % class
