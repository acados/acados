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

			% check for module name
			tokens = strsplit(field, '_');

			% symbolics
			if (strcmp(tokens{1}, 'sym'))

				if (strcmp(field, 'sym_x'))
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
				else
					disp(['acados_ocp_model: set: wrong field: ', field]);
				end

			% cost
			elseif (strcmp(tokens{1}, 'cost'))

				if (strcmp(field, 'cost_type'))
					obj.cost_type = value;
					obj.model_struct.cost_type = value;
				elseif (strcmp(field, 'cost_type_e'))
					obj.cost_e_type = value;
					obj.model_struct.cost_e_type = value;
				elseif (strcmp(field, 'cost_expr_y'))
					obj.expr_y = value;
					obj.model_struct.expr_y = value;
				elseif (strcmp(field, 'cost_param_y'))
					obj.param_y = value;
					obj.model_struct.param_y = value;
				elseif (strcmp(field, 'cost_expr_y_e'))
					obj.expr_y_e = value;
					obj.model_struct.expr_y_e = value;
				elseif (strcmp(field, 'cost_param_y_e'))
					obj.param_y_e = value;
					obj.model_struct.param_y_e = value;
				elseif (strcmp(field, 'cost_expr_ext_cost'))
					obj.expr_ext_cost = value;
					obj.model_struct.expr_ext_cost = value;
				elseif (strcmp(field, 'cost_param_ext_cost'))
					obj.param_ext_cost = value;
					obj.model_struct.param_ext_cost = value;
				elseif (strcmp(field, 'cost_expr_ext_cost_e'))
					obj.expr_ext_cost_e = value;
					obj.model_struct.expr_ext_cost_e = value;
				elseif (strcmp(field, 'cost_param_ext_cost_e'))
					obj.param_ext_cost_e = value;
					obj.model_struct.param_ext_cost_e = value;
				elseif (strcmp(field, 'cost_Vu'))
					obj.Vu = value;
					obj.model_struct.Vu = value;
				elseif (strcmp(field, 'cost_Vx'))
					obj.Vx = value;
					obj.model_struct.Vx = value;
				elseif (strcmp(field, 'cost_Vx_e'))
					obj.Vx_e = value;
					obj.model_struct.Vx_e = value;
				elseif (strcmp(field, 'cost_W'))
					obj.W = value;
					obj.model_struct.W = value;
				elseif (strcmp(field, 'cost_W_e'))
					obj.W_e = value;
					obj.model_struct.W_e = value;
				elseif (strcmp(field, 'cost_yr'))
					obj.yr = value;
					obj.model_struct.yr = value;
				elseif (strcmp(field, 'cost_yr_e'))
					obj.yr_e = value;
					obj.model_struct.yr_e = value;
				elseif (strcmp(field, 'cost_Z'))
					obj.Z = value;
					obj.model_struct.Z = value;
				elseif (strcmp(field, 'cost_Z_e'))
					obj.Z_e = value;
					obj.model_struct.Z_e = value;
				elseif (strcmp(field, 'cost_Zl'))
					obj.Zl = value;
					obj.model_struct.Zl = value;
				elseif (strcmp(field, 'cost_Zl_e'))
					obj.Zl_e = value;
					obj.model_struct.Zl_e = value;
				elseif (strcmp(field, 'cost_Zu'))
					obj.Zu = value;
					obj.model_struct.Zu = value;
				elseif (strcmp(field, 'cost_Zu_e'))
					obj.Zu_e = value;
					obj.model_struct.Zu_e = value;
				elseif (strcmp(field, 'cost_zl'))
					obj.zl = value;
					obj.model_struct.zl = value;
				elseif (strcmp(field, 'cost_zl_e'))
					obj.zl_e = value;
					obj.model_struct.zl_e = value;
				elseif (strcmp(field, 'cost_z'))
					obj.z = value;
					obj.model_struct.z = value;
				elseif (strcmp(field, 'cost_z_e'))
					obj.z_e = value;
					obj.model_struct.z_e = value;
				elseif (strcmp(field, 'cost_zu'))
					obj.zu = value;
					obj.model_struct.zu = value;
				elseif (strcmp(field, 'cost_zu_e'))
					obj.zu_e = value;
					obj.model_struct.zu_e = value;
				else
					disp(['acados_ocp_model: set: wrong field: ', field]);
				end

			% dynamics
			elseif (strcmp(tokens{1}, 'dyn'))

				if (strcmp(field, 'dyn_type'))
					obj.dyn_type = value;
					obj.model_struct.dyn_type = value;
				elseif (strcmp(field, 'dyn_expr_f'))
					obj.expr_f = value;
					obj.model_struct.expr_f = value;
				elseif (strcmp(field, 'dyn_param_f'))
					obj.param_f = value;
					obj.model_struct.param_f = value;
				else
					disp(['acados_ocp_model: set: wrong field: ', field]);
				end

			% constraints
			elseif (strcmp(tokens{1}, 'constr'))

				if (strcmp(field, 'constr_type'))
					obj.constr_type = value;
					obj.model_struct.constr_type = value;
				elseif (strcmp(field, 'constr_x0'))
					obj.x0 = value;
					obj.model_struct.x0 = value;
				elseif (strcmp(field, 'constr_Jbx'))
					obj.Jbx = value;
					obj.model_struct.Jbx = value;
				elseif (strcmp(field, 'constr_lbx'))
					obj.lbx = value;
					obj.model_struct.lbx = value;
				elseif (strcmp(field, 'constr_ubx'))
					obj.ubx = value;
					obj.model_struct.ubx = value;
				elseif (strcmp(field, 'constr_Jbu'))
					obj.Jbu = value;
					obj.model_struct.Jbu = value;
				elseif (strcmp(field, 'constr_lbu'))
					obj.lbu = value;
					obj.model_struct.lbu = value;
				elseif (strcmp(field, 'constr_ubu'))
					obj.ubu = value;
					obj.model_struct.ubu = value;
				elseif (strcmp(field, 'constr_C'))
					obj.C = value;
					obj.model_struct.C = value;
				elseif (strcmp(field, 'constr_D'))
					obj.D = value;
					obj.model_struct.D = value;
				elseif (strcmp(field, 'constr_lg'))
					obj.lg = value;
					obj.model_struct.lg = value;
				elseif (strcmp(field, 'constr_ug'))
					obj.ug = value;
					obj.model_struct.ug = value;
				elseif (strcmp(field, 'constr_C_e'))
					obj.C_e = value;
					obj.model_struct.C_e = value;
				elseif (strcmp(field, 'constr_lg_e'))
					obj.lg_e = value;
					obj.model_struct.lg_e = value;
				elseif (strcmp(field, 'constr_ug_e'))
					obj.ug_e = value;
					obj.model_struct.ug_e = value;
				elseif (strcmp(field, 'constr_expr_h'))
					obj.expr_h = value;
					obj.model_struct.expr_h = value;
				elseif (strcmp(field, 'constr_param_h'))
					obj.param_h = value;
					obj.model_struct.param_h = value;
				elseif (strcmp(field, 'constr_lh'))
					obj.lh = value;
					obj.model_struct.lh = value;
				elseif (strcmp(field, 'constr_uh'))
					obj.uh = value;
					obj.model_struct.uh = value;
				elseif (strcmp(field, 'constr_expr_h_e'))
					obj.expr_h_e = value;
					obj.model_struct.expr_h_e = value;
				elseif (strcmp(field, 'constr_param_h_e'))
					obj.param_h_e = value;
					obj.model_struct.param_h_e = value;
				elseif (strcmp(field, 'constr_lh_e'))
					obj.lh_e = value;
					obj.model_struct.lh_e = value;
				elseif (strcmp(field, 'constr_uh_e'))
					obj.uh_e = value;
					obj.model_struct.uh_e = value;
				elseif (strcmp(field, 'constr_Jsbu'))
					obj.Jsbu = value;
					obj.model_struct.Jsbu = value;
	%			elseif (strcmp(field, 'constr_lsbu'))
	%				obj.lsbu = value;
	%				obj.model_struct.lsbu = value;
	%			elseif (strcmp(field, 'constr_usbu'))
	%				obj.usbu = value;
	%				obj.model_struct.usbu = value;
				elseif (strcmp(field, 'constr_Jsbx'))
					obj.Jsbx = value;
					obj.model_struct.Jsbx = value;
	%			elseif (strcmp(field, 'constr_lsbx'))
	%				obj.lsbx = value;
	%				obj.model_struct.lsbx = value;
	%			elseif (strcmp(field, 'constr_usbx'))
	%				obj.usbx = value;
	%				obj.model_struct.usbx = value;
				elseif (strcmp(field, 'constr_Jsg'))
					obj.Jsg = value;
					obj.model_struct.Jsg = value;
	%			elseif (strcmp(field, 'constr_lsg'))
	%				obj.lsg = value;
	%				obj.model_struct.lsg = value;
	%			elseif (strcmp(field, 'constr_usg'))
	%				obj.usg = value;
	%				obj.model_struct.usg = value;
				elseif (strcmp(field, 'constr_Jsg_e'))
					obj.Jsg_e = value;
					obj.model_struct.Jsg_e = value;
	%			elseif (strcmp(field, 'constr_lsg_e'))
	%				obj.lsg_e = value;
	%				obj.model_struct.lsg_e = value;
	%			elseif (strcmp(field, 'constr_usg_e'))
	%				obj.usg_e = value;
	%				obj.model_struct.usg_e = value;
				elseif (strcmp(field, 'constr_Jsh'))
					obj.Jsh = value;
					obj.model_struct.Jsh = value;
	%			elseif (strcmp(field, 'constr_lsh'))
	%				obj.lsh = value;
	%				obj.model_struct.lsh = value;
	%			elseif (strcmp(field, 'constr_ush'))
	%				obj.ush = value;
	%				obj.model_struct.ush = value;
				elseif (strcmp(field, 'constr_Jsh_e'))
					obj.Jsh_e = value;
					obj.model_struct.Jsh_e = value;
	%			elseif (strcmp(field, 'constr_lsh_e'))
	%				obj.lsh_e = value;
	%				obj.model_struct.lsh_e = value;
	%			elseif (strcmp(field, 'constr_ush_e'))
	%				obj.ush_e = value;
	%				obj.model_struct.ush_e = value;
				else
					disp(['acados_ocp_model: set: wrong field: ', field]);
				end

			% dims
			elseif (strcmp(tokens{1}, 'dim'))

				if (strcmp(field, 'dim_nx'))
					obj.nx = value;
					obj.model_struct.nx = value;
				elseif (strcmp(field, 'dim_nu'))
					obj.nu = value;
					obj.model_struct.nu = value;
				elseif (strcmp(field, 'dim_nz'))
					obj.nz = value;
					obj.model_struct.nz = value;
				elseif (strcmp(field, 'dim_ny'))
					obj.ny = value;
					obj.model_struct.ny = value;
				elseif (strcmp(field, 'dim_ny_e'))
					obj.ny_e = value;
					obj.model_struct.ny_e = value;
				elseif (strcmp(field, 'dim_nbx'))
					obj.nbx = value;
					obj.model_struct.nbx = value;
				elseif (strcmp(field, 'dim_nbu'))
					obj.nbu = value;
					obj.model_struct.nbu = value;
				elseif (strcmp(field, 'dim_ng'))
					obj.ng = value;
					obj.model_struct.ng = value;
				elseif (strcmp(field, 'dim_ng_e'))
					obj.ng_e = value;
					obj.model_struct.ng_e = value;
				elseif (strcmp(field, 'dim_nh'))
					obj.nh = value;
					obj.model_struct.nh = value;
				elseif (strcmp(field, 'dim_nh_e'))
					obj.nh_e = value;
					obj.model_struct.nh_e = value;
				elseif (strcmp(field, 'dim_ns'))
					obj.ns = value;
					obj.model_struct.ns = value;
				elseif (strcmp(field, 'dim_ns_e'))
					obj.ns_e = value;
					obj.model_struct.ns_e = value;
				elseif (strcmp(field, 'dim_nsbu'))
					obj.nsbu = value;
					obj.model_struct.nsbu = value;
				elseif (strcmp(field, 'dim_nsbx'))
					obj.nsbx = value;
					obj.model_struct.nsbx = value;
				elseif (strcmp(field, 'dim_nsg'))
					obj.nsg = value;
					obj.model_struct.nsg = value;
				elseif (strcmp(field, 'dim_nsg_e'))
					obj.nsg_e = value;
					obj.model_struct.nsg_e = value;
				elseif (strcmp(field, 'dim_nsh'))
					obj.nsh = value;
					obj.model_struct.nsh = value;
				elseif (strcmp(field, 'dim_nsh_e'))
					obj.nsh_e = value;
					obj.model_struct.nsh_e = value;
				elseif (strcmp(field, 'dim_np'))
					obj.np = value;
					obj.model_struct.np = value;
				else
					disp(['acados_ocp_model: set: wrong field: ', field]);
				end

			% others
			else

				if (strcmp(field, 'T'))
					obj.T = value;
					obj.model_struct.T = value;

				else
					disp(['acados_ocp_model: set: wrong field: ', field]);
				end
			end
		end

	end % methods



end % class


