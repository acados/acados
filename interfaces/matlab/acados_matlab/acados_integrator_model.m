classdef acados_integrator_model < handle
	


	properties
		name
		type
		expr
		x
		u
		xdot
		z
		nx
		nu
		nz
		model_struct
	end %properties



	methods
		

		function obj = acados_integrator_model()
			obj.name = 'model';
			obj.type = 0;
			obj.expr = 0;
			obj.x = 0;
			obj.u = 0;
			obj.xdot = 0;
			obj.nx = 0;
			obj.nu = 0;
			obj.model_struct = struct;
			obj.model_struct.name = obj.name;
		end


		function obj = set(obj, field, value)
			if (strcmp(field, 'type'))
				obj.type = value;
				obj.model_struct.type = value;
			elseif (strcmp(field, 'expr'))
				obj.expr = value;
				obj.model_struct.expr = value;
			elseif (strcmp(field, 'x'))
				obj.x = value;
				obj.model_struct.x = value;
			elseif (strcmp(field, 'u'))
				obj.u = value;
				obj.model_struct.u = value;
			elseif (strcmp(field, 'xdot'))
				obj.xdot = value;
				obj.model_struct.xdot = value;
			elseif (strcmp(field, 'z'))
				obj.z = value;
				obj.model_struct.z = value;
			elseif (strcmp(field, 'nx'))
				obj.nx = value;
				obj.model_struct.nx = value;
			elseif (strcmp(field, 'nu'))
				obj.nu = value;
				obj.model_struct.nu = value;
			elseif (strcmp(field, 'nz'))
				obj.nz = value;
				obj.model_struct.nz = value;
			else
				disp('acados_integrator_model: set: wrong field');
			end
		end

	end % methods



end % class

