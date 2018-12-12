classdef acados_integrator_opts < handle
	


	properties
		codgen_model
		num_stages
		num_steps
		scheme
		sens_forw
		opts_struct
	end % properties



	methods
		

		function obj = acados_integrator_opts()
			obj.codgen_model = 'true';
			obj.num_stages = 0;
			obj.num_steps = 0;
			obj.scheme = 0;
			obj.sens_forw = 0;
			obj.opts_struct = struct;
			obj.opts_struct.codgen_model = obj.codgen_model;
		end


		function obj = set(obj, field, value)
			if (strcmp(field, 'codgen_model'))
				obj.codgen_model = value;
				obj.opts_struct.codgen_model = value;
			elseif (strcmp(field, 'num_stages'))
				obj.num_stages = value;
				obj.opts_struct.num_stages = value;
			elseif (strcmp(field, 'num_steps'))
				obj.num_steps = value;
				obj.opts_struct.num_steps = value;
			elseif (strcmp(field, 'scheme'))
				obj.scheme = value;
				obj.opts_struct.scheme = value;
			elseif (strcmp(field, 'sens_forw'))
				obj.sens_forw = value;
				if (strcmp(value, 'true'))
					obj.opts_struct.sens_forw = 1;
				else
					obj.opts_struct.sens_forw = 0;
				end
			else
				disp('acados_integrator_opts: set: wrong field');
			end
		end


	end % methods



end % class
