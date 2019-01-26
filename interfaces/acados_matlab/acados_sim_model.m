classdef acados_sim_model < handle
	


	properties
		name
		%dims
		nx
		nu
		nz
		np
		% symbolics
		sym_x
		sym_u
		sym_xdot
		sym_z
		sym_p
		% model
		dyn_type
		expr_f
		param_f
		T
		model_struct
	end %properties



	methods
		

		function obj = acados_sim_model()
			% default values
			obj.name = 'sim_model';
			% default values
			obj.model_struct = struct;
			obj.param_f = 'false';
			% initialize model struct
			obj.model_struct.name = obj.name;
			obj.model_struct.param_f = obj.param_f;
		end


		function obj = set(obj, field, value)
			% dims
			if (strcmp(field, 'nx'))
				obj.nx = value;
				obj.model_struct.nx = value;
			elseif (strcmp(field, 'nu'))
				obj.nu = value;
				obj.model_struct.nu = value;
			elseif (strcmp(field, 'nz'))
				obj.nz = value;
				obj.model_struct.nz = value;
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
			% model
			elseif (strcmp(field, 'dyn_type'))
				obj.dyn_type = value;
				obj.model_struct.dyn_type = value;
			elseif (strcmp(field, 'param_f'))
				obj.param_f = value;
				obj.model_struct.param_f = value;
			elseif (strcmp(field, 'expr_f'))
				obj.expr_f = value;
				obj.model_struct.expr_f = value;
			elseif (strcmp(field, 'T'))
				obj.T = value;
				obj.model_struct.T = value;
			else
				disp(['acados_sim_model: set: wrong field: ', field]);
			end
		end

	end % methods



end % class

