classdef acados_sim_model < handle
	


	properties
		name
		%dims
		dim_nx
		dim_nu
		dim_nz
		dim_np
		% symbolics
		sym_x
		sym_u
		sym_xdot
		sym_z
		sym_p
		% dynamics
		dyn_type
		dyn_expr_f
		dyn_param_f
		T
		model_struct
	end %properties



	methods
		

		function obj = acados_sim_model()
			% default values
			obj.name = 'sim_model';
			% default values
			obj.model_struct = struct;
			obj.dyn_param_f = 'false';
			% initialize model struct
			obj.model_struct.name = obj.name;
			obj.model_struct.dyn_param_f = obj.dyn_param_f;
		end


		function obj = set(obj, field, value)
			% dims
			if (strcmp(field, 'dim_nx'))
				obj.dim_nx = value;
				obj.model_struct.dim_nx = value;
			elseif (strcmp(field, 'dim_nu'))
				obj.dim_nu = value;
				obj.model_struct.dim_nu = value;
			elseif (strcmp(field, 'dim_nz'))
				obj.dim_nz = value;
				obj.model_struct.dim_nz = value;
			elseif (strcmp(field, 'dim_np'))
				obj.dim_np = value;
				obj.model_struct.dim_np = value;
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
			% dynamics
			elseif (strcmp(field, 'dyn_type'))
				obj.dyn_type = value;
				obj.model_struct.dyn_type = value;
			elseif (strcmp(field, 'dyn_param_f'))
				obj.dyn_param_f = value;
				obj.model_struct.dyn_param_f = value;
			elseif (strcmp(field, 'dyn_expr_f'))
				obj.dyn_expr_f = value;
				obj.model_struct.dyn_expr_f = value;
			elseif (strcmp(field, 'T'))
				obj.T = value;
				obj.model_struct.T = value;
			else
				disp(['acados_sim_model: set: wrong field: ', field]);
			end
		end

	end % methods



end % class

