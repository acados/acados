classdef acados_integrator_opts < handle
	
	properties
		acados
		py_opts
		scheme
		sens_forw
		codgen_model
	end


	methods
		
		function obj = acados_integrator_opts()
			obj.acados = py.importlib.import_module('acados');
			obj.py_opts = obj.acados.sim.acados_integrator_opts();
			obj.scheme = [];
			obj.sens_forw = 'false';
			obj.codgen_model = 'true';
			obj.py_opts.set('codgen_model', 'false');
		end

		function obj = set(obj, field, value)
			if (strcmp(field, 'scheme'))
				obj.scheme = value;
				obj.py_opts.set(field, value);
			elseif (strcmp(field, 'sens_forw'))
				obj.sens_forw = value;
				obj.py_opts.set(field, value);
			elseif (strcmp(field, 'codgen_model'))
				obj.codgen_model = value;
			end
		end
	
	end

end



