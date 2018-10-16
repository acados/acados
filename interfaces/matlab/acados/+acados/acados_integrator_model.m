classdef acados_integrator_model
	
	properties
		ac
		np
		py_acados_integrator_model
	end

	methods
		function obj = acados_integrator_model()
			obj.np = py.importlib.import_module('numpy');
			obj.ac = py.importlib.import_module('acados');
			obj.py_acados_integrator_model = obj.ac.sim.acados_integrator_model();
		end

		function set(obj, field, value)
            obj.py_acados_integrator_model.set(field, value);
		end
	end
end
