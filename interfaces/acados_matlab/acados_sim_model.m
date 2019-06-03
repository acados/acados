classdef acados_sim_model < handle
	


	properties
		model_struct
	end %properties



	methods
		

		function obj = acados_sim_model()
			obj.model_struct = struct;
			obj.model_struct.name = 'sim_model';
			obj.model_struct.dyn_param_f = 'false';
		end


		function obj = set(obj, field, value)
			% dims
			if (strcmp(field, 'dim_nx'))
				obj.model_struct.dim_nx = value;
			elseif (strcmp(field, 'dim_nu'))
				obj.model_struct.dim_nu = value;
			elseif (strcmp(field, 'dim_nz'))
				obj.model_struct.dim_nz = value;
			elseif (strcmp(field, 'dim_np'))
				obj.model_struct.dim_np = value;
			% symbolics
			elseif (strcmp(field, 'sym_x'))
				obj.model_struct.sym_x = value;
			elseif (strcmp(field, 'sym_xdot'))
				obj.model_struct.sym_xdot = value;
			elseif (strcmp(field, 'sym_u'))
				obj.model_struct.sym_u = value;
			elseif (strcmp(field, 'sym_z'))
				obj.model_struct.sym_z = value;
			elseif (strcmp(field, 'sym_p'))
				obj.model_struct.sym_p = value;
			% dynamics
			elseif (strcmp(field, 'dyn_type'))
				obj.model_struct.dyn_type = value;
			elseif (strcmp(field, 'dyn_param_f'))
				obj.model_struct.dyn_param_f = value;
			elseif (strcmp(field, 'dyn_expr_f'))
				obj.model_struct.dyn_expr_f = value;
			elseif (strcmp(field, 'T'))
				obj.model_struct.T = value;
			elseif (strcmp(field, 'seed_adj'))
				obj.model_struct.seed_adj = value;
			else
				disp(['acados_sim_model: set: wrong field: ', field]);
			end
		end

	end % methods



end % class

