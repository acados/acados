classdef acados_integrator_opts
	
	properties
		ac
		np
		py_acados_integrator_opts
	end

	methods
		function obj = acados_integrator_opts()
			obj.np = py.importlib.import_module('numpy');
			obj.ac = py.importlib.import_module('acados');
			obj.py_acados_integrator_opts = obj.ac.acados_integrator_opts();
		end

		function set(obj, field, value)
            % serialize casadi expressions
            if field == 'ode_expr'
                s_value = value.serialize()
                obj.py_acados_intergrator_model('serialized_set_in', 1)
            else
                s_value = value 
           end

            obj.py_acados_intergrator_opts.set(field, value);
		end
	end
end
