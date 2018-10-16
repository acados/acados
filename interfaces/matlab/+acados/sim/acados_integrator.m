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
			obj.py_acados_integrator_model = obj.ac.acados_integrator_model();
		end

		function set(obj, field, value)
            obj.py_acados_intergrator_model.set(field, value);
		end
	end
end

classdef acados_integrator
	
	properties
		ac
		np
		py_acados_integrator
	end

	methods
		function obj = acados_integrator(model, opts)
			obj.np = py.importlib.import_module('numpy');
			obj.ac = py.importlib.import_module('acados');
			obj.py_acados_integrator = obj.ac.acados_integrator(model.py_acados_model, opts.py_acados_otps);
		end

		function set(obj, field, value)
            obj.py_acados_intergrator.set(field, value);
		end

		function solve()
            obj.py_acados_intergrator.solve();
		end

		function get(obj, field, value)
            obj.py_acados_intergrator.get(field);
		end
	end
end
