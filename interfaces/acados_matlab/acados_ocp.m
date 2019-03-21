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

			if (strcmp(obj.opts_struct.codgen_model, 'true'))
				% select files to compile
				c_sources = ' ';
				% dynamics
				if (strcmp(obj.opts_struct.sim_method, 'erk'))
					% generate c for function and derivatives using casadi
					generate_c_code_explicit_ode(obj.model_struct, obj.opts_struct);
					% sources list
					c_sources = [c_sources, 'ocp_model_dyn_expl_ode_fun.c '];
					c_sources = [c_sources, 'ocp_model_dyn_expl_vde_for.c '];
					c_sources = [c_sources, 'ocp_model_dyn_expl_vde_adj.c '];
					c_sources = [c_sources, 'ocp_model_dyn_expl_ode_hes.c '];
				elseif (strcmp(obj.opts_struct.sim_method, 'irk'))
					% generate c for function and derivatives using casadi
					generate_c_code_implicit_ode(obj.model_struct, obj.opts_struct);
					% sources list
					c_sources = [c_sources, 'ocp_model_dyn_impl_ode_fun.c '];
					c_sources = [c_sources, 'ocp_model_dyn_impl_ode_fun_jac_x_xdot.c '];
					c_sources = [c_sources, 'ocp_model_dyn_impl_ode_fun_jac_x_xdot_u.c '];
					c_sources = [c_sources, 'ocp_model_dyn_impl_ode_jac_x_xdot_u.c '];
				else
					fprintf('\ncodegen_model: sim solver not supported: %s\n', obj.opts_struct.sim_method);
					return;
				end
				% nonlinear constraints
				if (strcmp(obj.model_struct.constr_type, 'bgh') && (isfield(obj.model_struct, 'constr_expr_h') || isfield(obj.model_struct, 'constr_expr_h_e')))
					% generate c for function and derivatives using casadi
					generate_c_code_nonlinear_constr(obj.model_struct, obj.opts_struct);
					% sources list
					if isfield(obj.model_struct, 'constr_expr_h')
						c_sources = [c_sources, 'ocp_model_constr_h_fun_jac_ut_xt.c '];
					end
					if isfield(obj.model_struct, 'constr_expr_h_e')
						c_sources = [c_sources, 'ocp_model_constr_h_e_fun_jac_ut_xt.c '];
					end
				end
				% nonlinear least squares
				if (strcmp(obj.model_struct.cost_type, 'nonlinear_ls') || strcmp(obj.model_struct.cost_type_e, 'nonlinear_ls'))
					% generate c for function and derivatives using casadi
					generate_c_code_nonlinear_least_squares(obj.model_struct, obj.opts_struct);
					% sources list
					if isfield(obj.model_struct, 'cost_expr_y')
						c_sources = [c_sources, 'ocp_model_cost_y_fun_jac_ut_xt.c '];
					end
					if isfield(obj.model_struct, 'cost_expr_y_e')
						c_sources = [c_sources, 'ocp_model_cost_y_e_fun_jac_ut_xt.c '];
					end
				end
				% external cost
				if (strcmp(obj.model_struct.cost_type, 'ext_cost') || strcmp(obj.model_struct.cost_type_e, 'ext_cost'))
					% generate c for function and derivatives using casadi
					generate_c_code_ext_cost(obj.model_struct, obj.opts_struct);
					% sources list
					if isfield(obj.model_struct, 'cost_expr_ext_cost')
						c_sources = [c_sources, 'ocp_model_cost_ext_cost_jac_hes.c '];
					end
					if isfield(obj.model_struct, 'cost_expr_ext_cost_e')
						c_sources = [c_sources, 'ocp_model_cost_ext_cost_e_jac_hes.c '];
					end
				end
				lib_name = ['libocp_model.so'];
				system(['gcc -O2 -fPIC -shared ', c_sources, ' -o ', lib_name]);
			end

			obj.C_ocp_ext_fun = ocp_create_ext_fun();

			% compile mex with model dependency
			if (strcmp(obj.opts_struct.compile_mex, 'true'))
				ocp_compile_mex_model_dep(obj.model_struct, obj.opts_struct);
			end

			% get pointers for external functions in model
			% dynamics
			if (strcmp(obj.opts_struct.sim_method, 'erk'))
				obj.C_ocp_ext_fun = ocp_set_ext_fun_dyn_expl(obj.C_ocp, obj.C_ocp_ext_fun, obj.model_struct, obj.opts_struct);
			elseif (strcmp(obj.opts_struct.sim_method, 'irk'))
				obj.C_ocp_ext_fun = ocp_set_ext_fun_dyn_impl(obj.C_ocp, obj.C_ocp_ext_fun, obj.model_struct, obj.opts_struct);
			else
				fprintf('\ncodegen_model: sim_method not supported: %s\n', obj.opts_struct.sim_method);
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

			% set in model
			ocp_set_model(obj.model_struct, obj.C_ocp_ext_fun, obj.C_ocp);
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

