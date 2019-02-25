classdef acados_ocp_model < handle
	


	properties
		name
		% dims
		T
		dim_nx
		dim_nu
		dim_nz
		dim_ny
		dim_ny_e
		dim_nbx
		dim_nbu
		dim_ng
		dim_ng_e
		dim_nh
		dim_nh_e
		dim_ns
		dim_ns_e
		dim_nsbu
		dim_nsbx
		dim_nsg
		dim_nsg_e
		dim_nsh
		dim_nsh_e
		dim_np
%		dim_np_e
		% symbolics
		sym_x
		sym_u
		sym_xdot
		sym_z
		sym_p
		% cost
		cost_type
		cost_type_e
		cost_expr_y
		cost_param_y
		cost_expr_y_e
		cost_param_y_e
		cost_expr_ext_cost
		cost_param_ext_cost
		cost_expr_ext_cost_e
		cost_param_ext_cost_e
		cost_Vu
		cost_Vx
		cost_Vx_e
		cost_W
		cost_W_e
		cost_yr
		cost_yr_e
		cost_Z
		cost_Z_e
		cost_Zl
		cost_Zl_e
		cost_Zu
		cost_Zu_e
		cost_z
		cost_z_e
		cost_zl
		cost_zl_e
		cost_zu
		cost_zu_e
		% dynamics
		dyn_type
		dyn_expr_f
		dyn_param_f
		% constraints
		constr_type
		constr_x0
		constr_Jbx
		constr_lbx
		constr_ubx
		constr_Jbu
		constr_lbu
		constr_ubu
		constr_C
		constr_D
		constr_lg
		constr_ug
		constr_C_e
		constr_lg_e
		constr_ug_e
		constr_expr_h
		constr_param_h
		constr_lh
		constr_uh
		constr_expr_h_e
		constr_param_h_e
		constr_lh_e
		constr_uh_e
		constr_Jsbu
%		constr_lsbu
%		constr_usbu
		constr_Jsbx
%		constr_lsbx
%		constr_usbx
		constr_Jsg
%		constr_lsg
%		constr_usg
		constr_Jsg_e
%		constr_lsg_e
%		constr_usg_e
		constr_Jsh
%		constr_lsh
%		constr_ush
		constr_Jsh_e
%		constr_lsh_e
%		constr_ush_e
		% structure
		model_struct
	end %properties



	methods
		

		function obj = acados_ocp_model()
			% default values
			obj.name = 'ocp_model';
			obj.cost_type = 'linear_ls';
			obj.cost_type_e = 'linear_ls';
			obj.cost_param_y = 'false';
			obj.cost_param_y_e = 'false';
			obj.cost_param_ext_cost = 'false';
			obj.cost_param_ext_cost_e = 'false';
			obj.dyn_type = 'implicit';
			obj.dyn_param_f = 'false';
			obj.constr_type = 'bgh';
			obj.constr_param_h = 'false';
			obj.constr_param_h_e = 'false';
			% model structure
			obj.model_struct = struct;
			% initialize model struct
			obj.model_struct.name = obj.name;
			obj.model_struct.cost_type = obj.cost_type;
			obj.model_struct.cost_type_e = obj.cost_type_e;
			obj.model_struct.cost_param_y = obj.cost_param_y;
			obj.model_struct.cost_param_y_e = obj.cost_param_y_e;
			obj.model_struct.cost_param_ext_cost = obj.cost_param_ext_cost;
			obj.model_struct.cost_param_ext_cost_e = obj.cost_param_ext_cost_e;
			obj.model_struct.dyn_type = obj.dyn_type;
			obj.model_struct.dyn_param_f = obj.dyn_param_f;
			obj.model_struct.constr_type = obj.constr_type;
			obj.model_struct.constr_param_h = obj.constr_param_h;
			obj.model_struct.constr_param_h_e = obj.constr_param_h_e;
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
					obj.cost_type_e = value;
					obj.model_struct.cost_type_e = value;
				elseif (strcmp(field, 'cost_expr_y'))
					obj.cost_expr_y = value;
					obj.model_struct.cost_expr_y = value;
				elseif (strcmp(field, 'cost_param_y'))
					obj.cost_param_y = value;
					obj.model_struct.cost_param_y = value;
				elseif (strcmp(field, 'cost_expr_y_e'))
					obj.cost_expr_y_e = value;
					obj.model_struct.cost_expr_y_e = value;
				elseif (strcmp(field, 'cost_param_y_e'))
					obj.cost_param_y_e = value;
					obj.model_struct.cost_param_y_e = value;
				elseif (strcmp(field, 'cost_expr_ext_cost'))
					obj.cost_expr_ext_cost = value;
					obj.model_struct.cost_expr_ext_cost = value;
				elseif (strcmp(field, 'cost_param_ext_cost'))
					obj.cost_param_ext_cost = value;
					obj.model_struct.cost_param_ext_cost = value;
				elseif (strcmp(field, 'cost_expr_ext_cost_e'))
					obj.cost_expr_ext_cost_e = value;
					obj.model_struct.cost_expr_ext_cost_e = value;
				elseif (strcmp(field, 'cost_param_ext_cost_e'))
					obj.cost_param_ext_cost_e = value;
					obj.model_struct.cost_param_ext_cost_e = value;
				elseif (strcmp(field, 'cost_Vu'))
					obj.cost_Vu = value;
					obj.model_struct.cost_Vu = value;
				elseif (strcmp(field, 'cost_Vx'))
					obj.cost_Vx = value;
					obj.model_struct.cost_Vx = value;
				elseif (strcmp(field, 'cost_Vx_e'))
					obj.cost_Vx_e = value;
					obj.model_struct.cost_Vx_e = value;
				elseif (strcmp(field, 'cost_W'))
					obj.cost_W = value;
					obj.model_struct.cost_W = value;
				elseif (strcmp(field, 'cost_W_e'))
					obj.cost_W_e = value;
					obj.model_struct.cost_W_e = value;
				elseif (strcmp(field, 'cost_yr'))
					obj.cost_yr = value;
					obj.model_struct.cost_yr = value;
				elseif (strcmp(field, 'cost_yr_e'))
					obj.cost_yr_e = value;
					obj.model_struct.cost_yr_e = value;
				elseif (strcmp(field, 'cost_Z'))
					obj.cost_Z = value;
					obj.model_struct.cost_Z = value;
				elseif (strcmp(field, 'cost_Z_e'))
					obj.cost_Z_e = value;
					obj.model_struct.cost_Z_e = value;
				elseif (strcmp(field, 'cost_Zl'))
					obj.cost_Zl = value;
					obj.model_struct.cost_Zl = value;
				elseif (strcmp(field, 'cost_Zl_e'))
					obj.cost_Zl_e = value;
					obj.model_struct.cost_Zl_e = value;
				elseif (strcmp(field, 'cost_Zu'))
					obj.cost_Zu = value;
					obj.model_struct.cost_Zu = value;
				elseif (strcmp(field, 'cost_Zu_e'))
					obj.cost_Zu_e = value;
					obj.model_struct.cost_Zu_e = value;
				elseif (strcmp(field, 'cost_zl'))
					obj.cost_zl = value;
					obj.model_struct.cost_zl = value;
				elseif (strcmp(field, 'cost_zl_e'))
					obj.cost_zl_e = value;
					obj.model_struct.cost_zl_e = value;
				elseif (strcmp(field, 'cost_z'))
					obj.cost_z = value;
					obj.model_struct.cost_z = value;
				elseif (strcmp(field, 'cost_z_e'))
					obj.cost_z_e = value;
					obj.model_struct.cost_z_e = value;
				elseif (strcmp(field, 'cost_zu'))
					obj.cost_zu = value;
					obj.model_struct.cost_zu = value;
				elseif (strcmp(field, 'cost_zu_e'))
					obj.cost_zu_e = value;
					obj.model_struct.cost_zu_e = value;
				else
					disp(['acados_ocp_model: set: wrong field: ', field]);
				end

			% dynamics
			elseif (strcmp(tokens{1}, 'dyn'))

				if (strcmp(field, 'dyn_type'))
					obj.dyn_type = value;
					obj.model_struct.dyn_type = value;
				elseif (strcmp(field, 'dyn_expr_f'))
					obj.dyn_expr_f = value;
					obj.model_struct.dyn_expr_f = value;
				elseif (strcmp(field, 'dyn_param_f'))
					obj.dyn_param_f = value;
					obj.model_struct.dyn_param_f = value;
				else
					disp(['acados_ocp_model: set: wrong field: ', field]);
				end

			% constraints
			elseif (strcmp(tokens{1}, 'constr'))

				if (strcmp(field, 'constr_type'))
					obj.constr_type = value;
					obj.model_struct.constr_type = value;
				elseif (strcmp(field, 'constr_x0'))
					obj.constr_x0 = value;
					obj.model_struct.constr_x0 = value;
				elseif (strcmp(field, 'constr_Jbx'))
					obj.constr_Jbx = value;
					obj.model_struct.constr_Jbx = value;
				elseif (strcmp(field, 'constr_lbx'))
					obj.constr_lbx = value;
					obj.model_struct.constr_lbx = value;
				elseif (strcmp(field, 'constr_ubx'))
					obj.constr_ubx = value;
					obj.model_struct.constr_ubx = value;
				elseif (strcmp(field, 'constr_Jbu'))
					obj.constr_Jbu = value;
					obj.model_struct.constr_Jbu = value;
				elseif (strcmp(field, 'constr_lbu'))
					obj.constr_lbu = value;
					obj.model_struct.constr_lbu = value;
				elseif (strcmp(field, 'constr_ubu'))
					obj.constr_ubu = value;
					obj.model_struct.constr_ubu = value;
				elseif (strcmp(field, 'constr_C'))
					obj.constr_C = value;
					obj.model_struct.constr_C = value;
				elseif (strcmp(field, 'constr_D'))
					obj.constr_D = value;
					obj.model_struct.constr_D = value;
				elseif (strcmp(field, 'constr_lg'))
					obj.constr_lg = value;
					obj.model_struct.constr_lg = value;
				elseif (strcmp(field, 'constr_ug'))
					obj.constr_ug = value;
					obj.model_struct.constr_ug = value;
				elseif (strcmp(field, 'constr_C_e'))
					obj.constr_C_e = value;
					obj.model_struct.constr_C_e = value;
				elseif (strcmp(field, 'constr_lg_e'))
					obj.constr_lg_e = value;
					obj.model_struct.constr_lg_e = value;
				elseif (strcmp(field, 'constr_ug_e'))
					obj.constr_ug_e = value;
					obj.model_struct.constr_ug_e = value;
				elseif (strcmp(field, 'constr_expr_h'))
					obj.constr_expr_h = value;
					obj.model_struct.constr_expr_h = value;
				elseif (strcmp(field, 'constr_param_h'))
					obj.constr_param_h = value;
					obj.model_struct.constr_param_h = value;
				elseif (strcmp(field, 'constr_lh'))
					obj.constr_lh = value;
					obj.model_struct.constr_lh = value;
				elseif (strcmp(field, 'constr_uh'))
					obj.constr_uh = value;
					obj.model_struct.constr_uh = value;
				elseif (strcmp(field, 'constr_expr_h_e'))
					obj.constr_expr_h_e = value;
					obj.model_struct.constr_expr_h_e = value;
				elseif (strcmp(field, 'constr_param_h_e'))
					obj.constr_param_h_e = value;
					obj.model_struct.constr_param_h_e = value;
				elseif (strcmp(field, 'constr_lh_e'))
					obj.constr_lh_e = value;
					obj.model_struct.constr_lh_e = value;
				elseif (strcmp(field, 'constr_uh_e'))
					obj.constr_uh_e = value;
					obj.model_struct.constr_uh_e = value;
				elseif (strcmp(field, 'constr_Jsbu'))
					obj.constr_Jsbu = value;
					obj.model_struct.constr_Jsbu = value;
	%			elseif (strcmp(field, 'constr_lsbu'))
	%				obj.constr_lsbu = value;
	%				obj.model_struct.constr_lsbu = value;
	%			elseif (strcmp(field, 'constr_usbu'))
	%				obj.constr_usbu = value;
	%				obj.model_struct.constr_usbu = value;
				elseif (strcmp(field, 'constr_Jsbx'))
					obj.constr_Jsbx = value;
					obj.model_struct.constr_Jsbx = value;
	%			elseif (strcmp(field, 'constr_lsbx'))
	%				obj.constr_lsbx = value;
	%				obj.model_struct.constr_lsbx = value;
	%			elseif (strcmp(field, 'constr_usbx'))
	%				obj.constr_usbx = value;
	%				obj.model_struct.constr_usbx = value;
				elseif (strcmp(field, 'constr_Jsg'))
					obj.constr_Jsg = value;
					obj.model_struct.constr_Jsg = value;
	%			elseif (strcmp(field, 'constr_lsg'))
	%				obj.constr_lsg = value;
	%				obj.model_struct.constr_lsg = value;
	%			elseif (strcmp(field, 'constr_usg'))
	%				obj.constr_usg = value;
	%				obj.model_struct.constr_usg = value;
				elseif (strcmp(field, 'constr_Jsg_e'))
					obj.constr_Jsg_e = value;
					obj.model_struct.constr_Jsg_e = value;
	%			elseif (strcmp(field, 'constr_lsg_e'))
	%				obj.constr_lsg_e = value;
	%				obj.model_struct.constr_lsg_e = value;
	%			elseif (strcmp(field, 'constr_usg_e'))
	%				obj.constr_usg_e = value;
	%				obj.model_struct.constr_usg_e = value;
				elseif (strcmp(field, 'constr_Jsh'))
					obj.constr_Jsh = value;
					obj.model_struct.constr_Jsh = value;
	%			elseif (strcmp(field, 'constr_lsh'))
	%				obj.constr_lsh = value;
	%				obj.model_struct.constr_lsh = value;
	%			elseif (strcmp(field, 'constr_ush'))
	%				obj.constr_ush = value;
	%				obj.model_struct.constr_ush = value;
				elseif (strcmp(field, 'constr_Jsh_e'))
					obj.constr_Jsh_e = value;
					obj.model_struct.constr_Jsh_e = value;
	%			elseif (strcmp(field, 'constr_lsh_e'))
	%				obj.constr_lsh_e = value;
	%				obj.model_struct.constr_lsh_e = value;
	%			elseif (strcmp(field, 'constr_ush_e'))
	%				obj.constr_ush_e = value;
	%				obj.model_struct.constr_ush_e = value;
				else
					disp(['acados_ocp_model: set: wrong field: ', field]);
				end

			% dims
			elseif (strcmp(tokens{1}, 'dim'))

				if (strcmp(field, 'dim_nx'))
					obj.dim_nx = value;
					obj.model_struct.dim_nx = value;
				elseif (strcmp(field, 'dim_nu'))
					obj.dim_nu = value;
					obj.model_struct.dim_nu = value;
				elseif (strcmp(field, 'dim_nz'))
					obj.dim_nz = value;
					obj.model_struct.dim_nz = value;
				elseif (strcmp(field, 'dim_ny'))
					obj.dim_ny = value;
					obj.model_struct.dim_ny = value;
				elseif (strcmp(field, 'dim_ny_e'))
					obj.dim_ny_e = value;
					obj.model_struct.dim_ny_e = value;
				elseif (strcmp(field, 'dim_nbx'))
					obj.dim_nbx = value;
					obj.model_struct.dim_nbx = value;
				elseif (strcmp(field, 'dim_nbu'))
					obj.dim_nbu = value;
					obj.model_struct.dim_nbu = value;
				elseif (strcmp(field, 'dim_ng'))
					obj.dim_ng = value;
					obj.model_struct.dim_ng = value;
				elseif (strcmp(field, 'dim_ng_e'))
					obj.dim_ng_e = value;
					obj.model_struct.dim_ng_e = value;
				elseif (strcmp(field, 'dim_nh'))
					obj.dim_nh = value;
					obj.model_struct.dim_nh = value;
				elseif (strcmp(field, 'dim_nh_e'))
					obj.dim_nh_e = value;
					obj.model_struct.dim_nh_e = value;
				elseif (strcmp(field, 'dim_ns'))
					obj.dim_ns = value;
					obj.model_struct.dim_ns = value;
				elseif (strcmp(field, 'dim_ns_e'))
					obj.dim_ns_e = value;
					obj.model_struct.dim_ns_e = value;
				elseif (strcmp(field, 'dim_nsbu'))
					obj.dim_nsbu = value;
					obj.model_struct.dim_nsbu = value;
				elseif (strcmp(field, 'dim_nsbx'))
					obj.dim_nsbx = value;
					obj.model_struct.dim_nsbx = value;
				elseif (strcmp(field, 'dim_nsg'))
					obj.dim_nsg = value;
					obj.model_struct.dim_nsg = value;
				elseif (strcmp(field, 'dim_nsg_e'))
					obj.dim_nsg_e = value;
					obj.model_struct.dim_nsg_e = value;
				elseif (strcmp(field, 'dim_nsh'))
					obj.dim_nsh = value;
					obj.model_struct.dim_nsh = value;
				elseif (strcmp(field, 'dim_nsh_e'))
					obj.dim_nsh_e = value;
					obj.model_struct.dim_nsh_e = value;
				elseif (strcmp(field, 'dim_np'))
					obj.dim_np = value;
					obj.model_struct.dim_np = value;
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


