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
		ng
		nh
		dyn_type
		dyn_expr
		constr_expr_h
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
		C
		D
		lg
		ug
		lh
		uh
		% structure
		model_struct
	end %properties



	methods
		

		function obj = acados_ocp_model()
			% default values
			obj.name = 'ocp_model';
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
			elseif (strcmp(field, 'ng'))
				obj.ng = value;
				obj.model_struct.ng = value;
			elseif (strcmp(field, 'nh'))
				obj.nh = value;
				obj.model_struct.nh = value;
			elseif (strcmp(field, 'dyn_type'))
				obj.dyn_type = value;
				obj.model_struct.dyn_type = value;
			elseif (strcmp(field, 'dyn_expr'))
				obj.dyn_expr = value;
				obj.model_struct.dyn_expr = value;
			elseif (strcmp(field, 'constr_expr_h'))
				obj.constr_expr_h = value;
				obj.model_struct.constr_expr_h = value;
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
			elseif (strcmp(field, 'C'))
				obj.C = value;
				obj.model_struct.C = value;
			elseif (strcmp(field, 'D'))
				obj.D = value;
				obj.model_struct.D = value;
			elseif (strcmp(field, 'lg'))
				obj.lg = value;
				obj.model_struct.lg = value;
			elseif (strcmp(field, 'ug'))
				obj.ug = value;
				obj.model_struct.ug = value;
			elseif (strcmp(field, 'lh'))
				obj.lh = value;
				obj.model_struct.lh = value;
			elseif (strcmp(field, 'uh'))
				obj.uh = value;
				obj.model_struct.uh = value;
			else
				disp(['acados_integrator_model: set: wrong field: ', field]);
			end
		end

	end % methods



end % class


