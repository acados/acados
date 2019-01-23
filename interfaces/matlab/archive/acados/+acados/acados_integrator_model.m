classdef acados_integrator_model < handle
	

	properties
		acados
		py_model
		model_name
		type
		ode_expr
		x
		u
		xdot
		z
		nx
		nu
		nz
		ode_expr_hash
	end


	methods
		
		function obj = acados_integrator_model()
%			obj.numpy = py.importlib.import_module('numpy');
			obj.acados = py.importlib.import_module('acados');
			obj.py_model = obj.acados.sim.acados_integrator_model();
			obj.model_name = []; % check with isempty()
			obj.type = [];
			obj.ode_expr = [];
			obj.x = [];
			obj.xdot = [];
			obj.u = [];
			obj.z = [];
			obj.nx = 0;
			obj.nu = 0;
			obj.nz = 0;
			obj.nz = [];
		end
		
		function obj = set(obj, field, value)
			if (strcmp(field, 'model_name'))
				obj.model_name = value;
				obj.py_model.set(field, value);
			elseif (strcmp(field, 'type'))
				obj.type = value;
				obj.py_model.set(field, value);
			elseif (strcmp(field, 'ode_expr'))
				obj.ode_expr = value;
				obj.ode_expr_hash = int64(py.hash(str(value)));
				obj.py_model.set('ode_expr_hash', int64(py.hash(str(value))));
			elseif (strcmp(field, 'ode_expr_hash'))
				obj.ode_expr_hash = value;
				obj.py_model.set('ode_expr_hash', value);
			elseif (strcmp(field, 'x'))
				obj.x = value;
				tmp = value.size();
				obj.nx = tmp(1);
				obj.py_model.set('nx', int32(obj.nx));
			elseif (strcmp(field, 'xdot'))
				obj.xdot = value;
			elseif (strcmp(field, 'u'))
				obj.u = value;
				tmp = value.size();
				obj.nu = tmp(1);
				obj.py_model.set('nu', int32(obj.nx));
			elseif (strcmp(field, 'z'))
				obj.z = value;
				tmp = value.size();
				obj.nz = tmp(1);
				obj.py_model.set('nz', int32(obj.nx));
			end

		end
	
	end

end


