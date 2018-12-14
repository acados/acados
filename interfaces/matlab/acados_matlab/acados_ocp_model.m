classdef acados_ocp_model < handle
	


	properties
		name
		T
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
		% cost
		Vul
		Vxl
		Vxm
		Wl
		Wm
		yrl
		yrm
		% constraints
		x0
		Jbx
		lbx
		ubx
		Jbu
		lbu
		ubu
		% structure
		model_struct
	end %properties



	methods
		

		function obj = acados_integrator_model()
			% default values
			obj.name = 'model';
			% model structure
			obj.model_struct = struct;
			% initialize model struct
			obj.model_struct.name = obj.name;
		end


		function obj = set(obj, field, value)
			if (strcmp(field, 'T'))
				obj.T = value;
				obj.model_struct.T = value;
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
			elseif (strcmp(field, 'Vul'))
				obj.Vul = value;
				obj.model_struct.Vul = value;
			elseif (strcmp(field, 'Vxl'))
				obj.Vxl = value;
				obj.model_struct.Vxl = value;
			elseif (strcmp(field, 'Vxm'))
				obj.Vxm = value;
				obj.model_struct.Vxm = value;
			elseif (strcmp(field, 'Wl'))
				obj.Wl = value;
				obj.model_struct.Wl = value;
			elseif (strcmp(field, 'Wm'))
				obj.Wm = value;
				obj.model_struct.Wm = value;
			elseif (strcmp(field, 'yrl'))
				obj.yrl = value;
				obj.model_struct.yrl = value;
			elseif (strcmp(field, 'yrm'))
				obj.yrm = value;
				obj.model_struct.yrm = value;
			elseif (strcmp(field, 'x0'))
				obj.x0 = value;
				obj.model_struct.x0 = value;
			elseif (strcmp(field, 'Jbx'))
				obj.Jbx = value;
				obj.model_struct.Jbx = value;
			elseif (strcmp(field, 'lbx'))
				obj.lbx = value;
				obj.model_struct.lbx = value;
			elseif (strcmp(field, 'ubx'))
				obj.ubx = value;
				obj.model_struct.ubx = value;
			elseif (strcmp(field, 'Jbu'))
				obj.Jbu = value;
				obj.model_struct.Jbu = value;
			elseif (strcmp(field, 'lbu'))
				obj.lbu = value;
				obj.model_struct.lbu = value;
			elseif (strcmp(field, 'ubu'))
				obj.ubu = value;
				obj.model_struct.ubu = value;
			else
				disp(['acados_integrator_model: set: wrong field: ', field]);
			end
		end

	end % methods



end % class


