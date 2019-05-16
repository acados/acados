classdef acados_ocp < handle

	properties
		C_ocp
		C_ocp_ext_fun
		model_struct
		opts_struct
	end % properties



	methods


		function obj = acados_ocp(model, opts)
			obj.model_struct = model.model_struct;
			obj.opts_struct = opts.opts_struct;

			% compile mex without model dependency
			if (strcmp(obj.opts_struct.compile_mex, 'true'))
				ocp_compile_mex();
			end

			obj.C_ocp = ocp_create(obj.model_struct, obj.opts_struct);

			% generate and compile casadi functions
			if (strcmp(obj.opts_struct.codgen_model, 'true'))
				ocp_compile_casadi_functions(obj.model_struct, obj.opts_struct);
			end

			obj.C_ocp_ext_fun = ocp_create_ext_fun();

			% compile mex with model dependency
			if (strcmp(obj.opts_struct.compile_mex, 'true'))
				ocp_compile_mex_model_dep(obj.model_struct, obj.opts_struct);
			end

			% get pointers for external functions in model
			% dynamics
			if (strcmp(obj.model_struct.dyn_type, 'explicit'))
				obj.C_ocp_ext_fun = ocp_set_ext_fun_dyn_expl(obj.C_ocp, obj.C_ocp_ext_fun, obj.model_struct, obj.opts_struct);
			elseif (strcmp(obj.model_struct.dyn_type, 'implicit'))
				obj.C_ocp_ext_fun = ocp_set_ext_fun_dyn_impl(obj.C_ocp, obj.C_ocp_ext_fun, obj.model_struct, obj.opts_struct);
			elseif (strcmp(obj.model_struct.dyn_type, 'discrete'))
				obj.C_ocp_ext_fun = ocp_set_ext_fun_dyn_disc(obj.C_ocp, obj.C_ocp_ext_fun, obj.model_struct, obj.opts_struct);
			else
				fprintf('\ncodegen_model: dyn_type not supported: %s\n', obj.model_struct.dyn_type);
				return;
			end
			% nonlinear constraints
			if (strcmp(obj.model_struct.constr_type, 'bgh') && isfield(obj.model_struct, 'constr_expr_h'))
				obj.C_ocp_ext_fun = ocp_set_ext_fun_constr_h(obj.C_ocp, obj.C_ocp_ext_fun, obj.model_struct, obj.opts_struct);
			end
			if (strcmp(obj.model_struct.constr_type, 'bgh') && isfield(obj.model_struct, 'constr_expr_h_e'))
				obj.C_ocp_ext_fun = ocp_set_ext_fun_constr_h_e(obj.C_ocp, obj.C_ocp_ext_fun, obj.model_struct, obj.opts_struct);
			end
			% nonlinear least squares
			if (strcmp(obj.model_struct.cost_type, 'nonlinear_ls'))
				obj.C_ocp_ext_fun = ocp_set_ext_fun_cost_y(obj.C_ocp, obj.C_ocp_ext_fun, obj.model_struct, obj.opts_struct);
			end
			if (strcmp(obj.model_struct.cost_type_e, 'nonlinear_ls'))
				obj.C_ocp_ext_fun = ocp_set_ext_fun_cost_y_e(obj.C_ocp, obj.C_ocp_ext_fun, obj.model_struct, obj.opts_struct);
			end
			% external cost
			if (strcmp(obj.model_struct.cost_type, 'ext_cost'))
				obj.C_ocp_ext_fun = ocp_set_ext_fun_cost_ext_cost(obj.C_ocp, obj.C_ocp_ext_fun, obj.model_struct, obj.opts_struct);
			end
			if (strcmp(obj.model_struct.cost_type_e, 'ext_cost'))
				obj.C_ocp_ext_fun = ocp_set_ext_fun_cost_ext_cost_e(obj.C_ocp, obj.C_ocp_ext_fun, obj.model_struct, obj.opts_struct);
			end
		end



		function solve(obj)
			ocp_solve(obj.C_ocp);
		end



%		function set(obj, field, value)
%			ocp_set(obj.model_struct, obj.opts_struct, obj.C_ocp, obj.C_ocp_ext_fun, field, value);
%		end
		function set(varargin)
			if nargin==3
				obj = varargin{1};
				field = varargin{2};
				value = varargin{3};
				ocp_set(obj.model_struct, obj.opts_struct, obj.C_ocp, obj.C_ocp_ext_fun, field, value);
			elseif nargin==4
				obj = varargin{1};
				field = varargin{2};
				stage = varargin{3};
				value = varargin{4};
				ocp_set(obj.model_struct, obj.opts_struct, obj.C_ocp, obj.C_ocp_ext_fun, field, value, stage);
			else
				disp('acados_ocp.set: wrong number of input arguments (2 or 3 allowed)');
			end
		end



		function value = get(varargin)
			if nargin==2
				obj = varargin{1};
				field = varargin{2};
				value = ocp_get(obj.C_ocp, field);
			elseif nargin==3
				obj = varargin{1};
				field = varargin{2};
				stage = varargin{3};
				value = ocp_get(obj.C_ocp, field, stage);
			else
				disp('acados_ocp.get: wrong number of input arguments (1 or 2 allowed)');
			end
		end



		function delete(obj)
			ocp_destroy_ext_fun(obj.model_struct, obj.C_ocp, obj.C_ocp_ext_fun);
			ocp_destroy(obj.C_ocp);
		end


	end % methods



end % class

