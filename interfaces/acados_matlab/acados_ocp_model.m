classdef acados_ocp_model < handle
	


	properties
		name
		% dims
		T
		nx
		nu
		nz
		ny
		ny_e
		nbx
		nbu
		ng
		ng_e
		nh
		nh_e
		ns
		ns_e
		nsbu
		nsbx
		nsg
		nsg_e
		nsh
		nsh_e
		np
%		np_e
		% symbolics
		sym_x
		sym_u
		sym_xdot
		sym_z
		sym_p
		% cost
		cost_type
		cost_e_type
		expr_y
		param_y
		expr_y_e
		param_y_e
		expr_ext_cost
		param_ext_cost
		expr_ext_cost_e
		param_ext_cost_e
		Vu
		Vx
		Vx_e
		W
		W_e
		yr
		yr_e
		Z
		Z_e
		Zl
		Zl_e
		Zu
		Zu_e
		z
		z_e
		zl
		zl_e
		zu
		zu_e
		% dynamics
		dyn_type
		expr_f
		param_f
		% constraints
		constr_type
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
		C_e
		lg_e
		ug_e
		expr_h
		param_h
		lh
		uh
		expr_h_e
		param_h_e
		lh_e
		uh_e
		Jsbu
%		lsbu
%		usbu
		Jsbx
%		lsbx
%		usbx
		Jsg
%		lsg
%		usg
		Jsg_e
%		lsg_e
%		usg_e
		Jsh
%		lsh
%		ush
		Jsh_e
%		lsh_e
%		ush_e
		% structure
		model_struct
	end %properties



	methods
		

		function obj = acados_ocp_model()
			% default values
			obj.name = 'ocp_model';
			obj.cost_type = 'linear_ls';
			obj.cost_e_type = 'linear_ls';
			obj.dyn_type = 'implicit';
			obj.constr_type = 'bgh';
			obj.param_y = 'false';
			obj.param_y_e = 'false';
			obj.param_f = 'false';
			obj.param_h = 'false';
			obj.param_h_e = 'false';
			obj.param_ext_cost = 'false';
			obj.param_ext_cost_e = 'false';
			% model structure
			obj.model_struct = struct;
			% initialize model struct
			obj.model_struct.name = obj.name;
			obj.model_struct.cost_type = obj.cost_type;
			obj.model_struct.cost_e_type = obj.cost_type;
			obj.model_struct.dyn_type = obj.dyn_type;
			obj.model_struct.constr_type = obj.constr_type;
			obj.model_struct.param_y = obj.param_y;
			obj.model_struct.param_y_e = obj.param_y_e;
			obj.model_struct.param_f = obj.param_f;
			obj.model_struct.param_h = obj.param_h;
			obj.model_struct.param_h_e = obj.param_h_e;
			obj.model_struct.param_ext_cost = obj.param_ext_cost;
			obj.model_struct.param_ext_cost_e = obj.param_ext_cost_e;
		end


		function obj = set(obj, field, value)
			% dims
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
			elseif (strcmp(field, 'ny'))
				obj.ny = value;
				obj.model_struct.ny = value;
			elseif (strcmp(field, 'ny_e'))
				obj.ny_e = value;
				obj.model_struct.ny_e = value;
			elseif (strcmp(field, 'nbx'))
				obj.nbx = value;
				obj.model_struct.nbx = value;
			elseif (strcmp(field, 'nbu'))
				obj.nbu = value;
				obj.model_struct.nbu = value;
			elseif (strcmp(field, 'ng'))
				obj.ng = value;
				obj.model_struct.ng = value;
			elseif (strcmp(field, 'ng_e'))
				obj.ng_e = value;
				obj.model_struct.ng_e = value;
			elseif (strcmp(field, 'nh'))
				obj.nh = value;
				obj.model_struct.nh = value;
			elseif (strcmp(field, 'nh_e'))
				obj.nh_e = value;
				obj.model_struct.nh_e = value;
			elseif (strcmp(field, 'ns'))
				obj.ns = value;
				obj.model_struct.ns = value;
			elseif (strcmp(field, 'ns_e'))
				obj.ns_e = value;
				obj.model_struct.ns_e = value;
			elseif (strcmp(field, 'nsbu'))
				obj.nsbu = value;
				obj.model_struct.nsbu = value;
			elseif (strcmp(field, 'nsbx'))
				obj.nsbx = value;
				obj.model_struct.nsbx = value;
			elseif (strcmp(field, 'nsg'))
				obj.nsg = value;
				obj.model_struct.nsg = value;
			elseif (strcmp(field, 'nsg_e'))
				obj.nsg_e = value;
				obj.model_struct.nsg_e = value;
			elseif (strcmp(field, 'nsh'))
				obj.nsh = value;
				obj.model_struct.nsh = value;
			elseif (strcmp(field, 'nsh_e'))
				obj.nsh_e = value;
				obj.model_struct.nsh_e = value;
			elseif (strcmp(field, 'np'))
				obj.np = value;
				obj.model_struct.np = value;
			% symbolics
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
			elseif (strcmp(field, 'sym_p'))
				obj.sym_p = value;
				obj.model_struct.sym_p = value;
			% cost
			elseif (strcmp(field, 'cost_type'))
				obj.cost_type = value;
				obj.model_struct.cost_type = value;
			elseif (strcmp(field, 'cost_e_type'))
				obj.cost_e_type = value;
				obj.model_struct.cost_e_type = value;
			elseif (strcmp(field, 'expr_y'))
				obj.expr_y = value;
				obj.model_struct.expr_y = value;
			elseif (strcmp(field, 'param_y'))
				obj.param_y = value;
				obj.model_struct.param_y = value;
			elseif (strcmp(field, 'expr_y_e'))
				obj.expr_y_e = value;
				obj.model_struct.expr_y_e = value;
			elseif (strcmp(field, 'param_y_e'))
				obj.param_y_e = value;
				obj.model_struct.param_y_e = value;
			elseif (strcmp(field, 'expr_ext_cost'))
				obj.expr_ext_cost = value;
				obj.model_struct.expr_ext_cost = value;
			elseif (strcmp(field, 'param_ext_cost'))
				obj.param_ext_cost = value;
				obj.model_struct.param_ext_cost = value;
			elseif (strcmp(field, 'expr_ext_cost_e'))
				obj.expr_ext_cost_e = value;
				obj.model_struct.expr_ext_cost_e = value;
			elseif (strcmp(field, 'param_ext_cost_e'))
				obj.param_ext_cost_e = value;
				obj.model_struct.param_ext_cost_e = value;
			elseif (strcmp(field, 'Vu'))
				obj.Vu = value;
				obj.model_struct.Vu = value;
			elseif (strcmp(field, 'Vx'))
				obj.Vx = value;
				obj.model_struct.Vx = value;
			elseif (strcmp(field, 'Vx_e'))
				obj.Vx_e = value;
				obj.model_struct.Vx_e = value;
			elseif (strcmp(field, 'W'))
				obj.W = value;
				obj.model_struct.W = value;
			elseif (strcmp(field, 'W_e'))
				obj.W_e = value;
				obj.model_struct.W_e = value;
			elseif (strcmp(field, 'yr'))
				obj.yr = value;
				obj.model_struct.yr = value;
			elseif (strcmp(field, 'yr_e'))
				obj.yr_e = value;
				obj.model_struct.yr_e = value;
			elseif (strcmp(field, 'Z'))
				obj.Z = value;
				obj.model_struct.Z = value;
			elseif (strcmp(field, 'Z_e'))
				obj.Z_e = value;
				obj.model_struct.Z_e = value;
			elseif (strcmp(field, 'Zl'))
				obj.Zl = value;
				obj.model_struct.Zl = value;
			elseif (strcmp(field, 'Zl_e'))
				obj.Zl_e = value;
				obj.model_struct.Zl_e = value;
			elseif (strcmp(field, 'Zu'))
				obj.Zu = value;
				obj.model_struct.Zu = value;
			elseif (strcmp(field, 'Zu_e'))
				obj.Zu_e = value;
				obj.model_struct.Zu_e = value;
			elseif (strcmp(field, 'zl'))
				obj.zl = value;
				obj.model_struct.zl = value;
			elseif (strcmp(field, 'zl_e'))
				obj.zl_e = value;
				obj.model_struct.zl_e = value;
			elseif (strcmp(field, 'z'))
				obj.z = value;
				obj.model_struct.z = value;
			elseif (strcmp(field, 'z_e'))
				obj.z_e = value;
				obj.model_struct.z_e = value;
			elseif (strcmp(field, 'zu'))
				obj.zu = value;
				obj.model_struct.zu = value;
			elseif (strcmp(field, 'zu_e'))
				obj.zu_e = value;
				obj.model_struct.zu_e = value;
			% dynamics
			elseif (strcmp(field, 'dyn_type'))
				obj.dyn_type = value;
				obj.model_struct.dyn_type = value;
			elseif (strcmp(field, 'expr_f'))
				obj.expr_f = value;
				obj.model_struct.expr_f = value;
			elseif (strcmp(field, 'param_f'))
				obj.param_f = value;
				obj.model_struct.param_f = value;
			% constraints
			elseif (strcmp(field, 'constr_type'))
				obj.constr_type = value;
				obj.model_struct.constr_type = value;
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
			elseif (strcmp(field, 'C_e'))
				obj.C_e = value;
				obj.model_struct.C_e = value;
			elseif (strcmp(field, 'lg_e'))
				obj.lg_e = value;
				obj.model_struct.lg_e = value;
			elseif (strcmp(field, 'ug_e'))
				obj.ug_e = value;
				obj.model_struct.ug_e = value;
			elseif (strcmp(field, 'expr_h'))
				obj.expr_h = value;
				obj.model_struct.expr_h = value;
			elseif (strcmp(field, 'param_h'))
				obj.param_h = value;
				obj.model_struct.param_h = value;
			elseif (strcmp(field, 'lh'))
				obj.lh = value;
				obj.model_struct.lh = value;
			elseif (strcmp(field, 'uh'))
				obj.uh = value;
				obj.model_struct.uh = value;
			elseif (strcmp(field, 'expr_h_e'))
				obj.expr_h_e = value;
				obj.model_struct.expr_h_e = value;
			elseif (strcmp(field, 'param_h_e'))
				obj.param_h_e = value;
				obj.model_struct.param_h_e = value;
			elseif (strcmp(field, 'lh_e'))
				obj.lh_e = value;
				obj.model_struct.lh_e = value;
			elseif (strcmp(field, 'uh_e'))
				obj.uh_e = value;
				obj.model_struct.uh_e = value;
			elseif (strcmp(field, 'Jsbu'))
				obj.Jsbu = value;
				obj.model_struct.Jsbu = value;
%			elseif (strcmp(field, 'lsbu'))
%				obj.lsbu = value;
%				obj.model_struct.lsbu = value;
%			elseif (strcmp(field, 'usbu'))
%				obj.usbu = value;
%				obj.model_struct.usbu = value;
			elseif (strcmp(field, 'Jsbx'))
				obj.Jsbx = value;
				obj.model_struct.Jsbx = value;
%			elseif (strcmp(field, 'lsbx'))
%				obj.lsbx = value;
%				obj.model_struct.lsbx = value;
%			elseif (strcmp(field, 'usbx'))
%				obj.usbx = value;
%				obj.model_struct.usbx = value;
			elseif (strcmp(field, 'Jsg'))
				obj.Jsg = value;
				obj.model_struct.Jsg = value;
%			elseif (strcmp(field, 'lsg'))
%				obj.lsg = value;
%				obj.model_struct.lsg = value;
%			elseif (strcmp(field, 'usg'))
%				obj.usg = value;
%				obj.model_struct.usg = value;
			elseif (strcmp(field, 'Jsg_e'))
				obj.Jsg_e = value;
				obj.model_struct.Jsg_e = value;
%			elseif (strcmp(field, 'lsg_e'))
%				obj.lsg_e = value;
%				obj.model_struct.lsg_e = value;
%			elseif (strcmp(field, 'usg_e'))
%				obj.usg_e = value;
%				obj.model_struct.usg_e = value;
			elseif (strcmp(field, 'Jsh'))
				obj.Jsh = value;
				obj.model_struct.Jsh = value;
%			elseif (strcmp(field, 'lsh'))
%				obj.lsh = value;
%				obj.model_struct.lsh = value;
%			elseif (strcmp(field, 'ush'))
%				obj.ush = value;
%				obj.model_struct.ush = value;
			elseif (strcmp(field, 'Jsh_e'))
				obj.Jsh_e = value;
				obj.model_struct.Jsh_e = value;
%			elseif (strcmp(field, 'lsh_e'))
%				obj.lsh_e = value;
%				obj.model_struct.lsh_e = value;
%			elseif (strcmp(field, 'ush_e'))
%				obj.ush_e = value;
%				obj.model_struct.ush_e = value;
			else
				disp(['acados_integrator_model: set: wrong field: ', field]);
			end
		end

	end % methods



end % class


