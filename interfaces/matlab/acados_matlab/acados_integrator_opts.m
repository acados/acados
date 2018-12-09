classdef acados_integrator_opts < handle
	


	properties
		codgen_model
		num_stages
		num_steps
		scheme
		sens_forw
	end % properties



	methods
		

		function obj = acados_integrator_opts()
			obj.codgen_model = 'true';
			obj.num_stages = 4;
			obj.num_steps = 5;
			obj.scheme = 'erk';
			obj.sens_forw = 'false';
		end


		function obj = set(obj, field, value)
			if (strcmp(field, 'codgen_model'))
				obj.codgen_model = value;
			elseif (strcmp(field, 'num_stages'))
				obj.num_stages = value;
			elseif (strcmp(field, 'num_steps'))
				obj.num_steps = value;
			elseif (strcmp(field, 'scheme'))
				obj.scheme = value;
			elseif (strcmp(field, 'sens_forw'))
				obj.sens_forw = value;
			end
		end


	end % methods



end % class
