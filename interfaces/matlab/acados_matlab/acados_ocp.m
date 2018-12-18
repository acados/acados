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
					c_sources = [c_sources, 'ocp_model_expl_ode_fun.c '];
					c_sources = [c_sources, 'ocp_model_expl_vde_for.c '];
					c_sources = [c_sources, 'ocp_model_expl_vde_adj.c '];
					c_sources = [c_sources, 'ocp_model_expl_ode_hes.c '];
				elseif (strcmp(obj.opts_struct.sim_method, 'irk'))
					% generate c for function and derivatives using casadi
					generate_c_code_implicit_ode(obj.model_struct, obj.opts_struct);
					% sources list
					c_sources = [c_sources, 'ocp_model_impl_ode_fun.c '];
					c_sources = [c_sources, 'ocp_model_impl_ode_fun_jac_x_xdot.c '];
					c_sources = [c_sources, 'ocp_model_impl_ode_fun_jac_x_xdot_u.c '];
					c_sources = [c_sources, 'ocp_model_impl_ode_jac_x_xdot_u.c '];
				else
					fprintf('\ncodegen_model: sim solver not supported: %s\n', obj.opts_struct.sim_method);
					return;
				end
				% nonlinear constraints
				if ((isfield(obj.model_struct, 'nh') && obj.model_struct.nh>0))
					% generate c for function and derivatives using casadi
					generate_c_code_nonlinear_constr(obj.model_struct, obj.opts_struct);
					% sources list
					c_sources = [c_sources, 'ocp_model_h_fun_jac_ut_xt.c '];
				end
				lib_name = ['libocp_model.so'];
				system(['gcc -fPIC -shared ', c_sources, ' -o ', lib_name]);
			end

			obj.C_ocp_ext_fun = ocp_create_ext_fun();

			% compile mex with model dependency
			if (strcmp(obj.opts_struct.compile_mex, 'true'))
				ocp_compile_mex_model_dep(obj.model_struct, obj.opts_struct);
			end

			% get pointers for external functions in model
			% dynamics
			if (strcmp(obj.opts_struct.sim_method, 'erk'))
				ocp_set_ext_fun_expl(obj.C_ocp_ext_fun, obj.model_struct, obj.opts_struct);
			elseif (strcmp(obj.opts_struct.sim_method, 'irk'))
				ocp_set_ext_fun_impl(obj.C_ocp_ext_fun, obj.model_struct, obj.opts_struct);
			else
				fprintf('\ncodegen_model: sim_method not supported: %s\n', obj.opts_struct.sim_method);
				return;
			end
			% nonlinear constraints
			if ((isfield(obj.model_struct, 'nh') && obj.model_struct.nh>0))
				ocp_set_ext_fun_h(obj.C_ocp_ext_fun, obj.model_struct, obj.opts_struct);
			end

			% set in model
			ocp_set_model(obj.C_ocp_ext_fun, obj.C_ocp);
		end



		function solve(obj)
			ocp_solve(obj.C_ocp);
		end



		function set(obj, field, value)
			ocp_set(obj.C_ocp, field, value);
		end



		function value = get(obj, field)
			value = ocp_get(obj.C_ocp, field);
		end



		function delete(obj)
%			fprintf('\nin delete\n');
			ocp_destroy(obj.C_ocp);
			ocp_destroy_ext_fun(obj.C_ocp_ext_fun);
		end


	end % methods



end % class

