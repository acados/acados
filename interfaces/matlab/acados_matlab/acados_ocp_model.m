classdef acados_ocp_model < handle
	


	properties
		name
		param_type
		N
		Ts
		nx
		nu
		nz
		nyl
		nym
		nbx
		nbu
		dyn_type
		dyn_expr
		sym_x
		sym_u
		sym_xdot
		sym_z
		Wl
		Wm
		model_struct
	end %properties



	methods
		

		function obj = acados_integrator_model()
			obj.name = 'model';
			obj.param_type = 'mult_shoot';
			obj.N = 0;
			obj.Ts = 0; % (uniform) sampling time
			obj.nx = 0;
			obj.nu = 0;
			obj.nz = 0;
			obj.nyl = 0;
			obj.nym = 0;
			obj.nbx = 0;
			obj.nux = 0;
			obj.dyn_type = 0;
			obj.dyn_expr = 0;
			obj.sym_x = 0;
			obj.sym_u = 0;
			obj.sym_xdot = 0;
			obj.sym_z = 0;
			obj.Wl = 0;
			obj.Wm = 0;
			% model structure
			obj.model_struct = struct;
			% fixed field values
			obj.model_struct.name = obj.name;
			obj.model_struct.param_type = obj.param_type;
		end


		function obj = set(obj, field, value)
			if (strcmp(field, 'N'))
				obj.N = value;
				obj.model_struct.N = value;
			elseif (strcmp(field, 'Ts'))
				obj.Ts = value;
				obj.model_struct.Ts = value;
			elseif (strcmp(field, 'nx'))
				obj.nx = value;
				obj.model_struct.nx = value;
			elseif (strcmp(field, 'nu'))
				obj.nu = value;
				obj.model_struct.nu = value;
			elseif (strcmp(field, 'nz'))
				obj.nz = value;
				obj.model_struct.nz = value;
			elseif (strcmp(field, 'nyl'))
				obj.nyl = value;
				obj.model_struct.nyl = value;
			elseif (strcmp(field, 'nym'))
				obj.nym = value;
				obj.model_struct.nym = value;
			elseif (strcmp(field, 'nbx'))
				obj.nbx = value;
				obj.model_struct.nbx = value;
			elseif (strcmp(field, 'nbu'))
				obj.nbu = value;
				obj.model_struct.nbu = value;
			elseif (strcmp(field, 'T'))
				obj.T = value;
				obj.model_struct.T = value;
			elseif (strcmp(field, 'dyn_type'))
				obj.type = value;
				obj.model_struct.type = value;
			elseif (strcmp(field, 'dyn_expr'))
				obj.expr = value;
				obj.model_struct.expr = value;
			elseif (strcmp(field, 'sym_x'))
				obj.sym_x = value;
				obj.model_struct.sym_x = value;
			elseif (strcmp(field, 'sym_xdot'))
				obj.sym_xdot = value;
				obj.model_struct.sym_xdot = value;
			elseif (strcmp(field, 'sym_u'))
				obj.sym_u = value;
				obj.model_struct.sym_u = value;
			elseif (strcmp(field, 'sym_z'))
				obj.sym_z = value;
				obj.model_struct.sym_z = value;
			elseif (strcmp(field, 'Wl'))
				obj.Wl = value;
				obj.model_struct.Wl = value;
			elseif (strcmp(field, 'Wm'))
				obj.Wm = value;
				obj.model_struct.Wm = value;
			else
				disp('acados_integrator_model: set: wrong field');
			end
		end

	end % methods



end % class


