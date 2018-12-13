classdef acados_ocp_opts < handle
	


	properties
		codgen_model
		sim_scheme
		opts_struct
	end % properties



	methods
		

		function obj = acados_ocp_opts()
			obj.codgen_model = 'true';
			obj.sim_scheme = 0;
			obj.opts_struct = struct;
			obj.opts_struct.codgen_model = obj.codgen_model;
		end


		function obj = set(obj, field, value)
			if (strcmp(field, 'codgen_model'))
				obj.codgen_model = value;
				obj.opts_struct.codgen_model = value;
			elseif (strcmp(field, 'sim_scheme'))
				obj.sim_scheme = value;
				obj.opts_struct.sim_scheme = value;
			else
				disp('acados_ocp_opts: set: wrong field');
			end
		end


	end % methods



end % class
