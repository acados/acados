classdef acados_sim < handle

	properties
		C_sim
		C_sim_ext_fun
		model_struct
		opts_struct
	end % properties



	methods


		function obj = acados_sim(model, opts)
			obj.model_struct = model.model_struct;
			obj.opts_struct = opts.opts_struct;

			% compile mex without model dependency
			if (strcmp(obj.opts_struct.compile_mex, 'true'))
				sim_compile_mex();
			end

			obj.C_sim = sim_create(obj.model_struct, obj.opts_struct);

			if obj.opts_struct.codgen_model
				c_sources = ' ';
				if (strcmp(obj.opts_struct.method, 'erk'))
					% generate c for function and derivatives using casadi
					generate_c_code_explicit_ode(obj.model_struct);
					% compile the code in a shared library
					c_sources = [c_sources, 'sim_model_dyn_expl_ode_fun.c '];
					c_sources = [c_sources, 'sim_model_dyn_expl_vde_for.c '];
					c_sources = [c_sources, 'sim_model_dyn_expl_vde_adj.c '];
					c_sources = [c_sources, 'sim_model_dyn_expl_ode_hes.c '];
				elseif (strcmp(obj.opts_struct.method, 'irk'))
					% generate c for function and derivatives using casadi
					generate_c_code_implicit_ode(obj.model_struct);
					% compile the code in a shared library
					c_sources = [c_sources, 'sim_model_dyn_impl_ode_fun.c '];
					c_sources = [c_sources, 'sim_model_dyn_impl_ode_fun_jac_x_xdot.c '];
					c_sources = [c_sources, 'sim_model_dyn_impl_ode_fun_jac_x_xdot_u.c '];
					c_sources = [c_sources, 'sim_model_dyn_impl_ode_jac_x_xdot_u.c '];
				else
					fprintf('\ncodegen_model: method not supported: %s\n', obj.opts_struct.method);
					return;
				end
				lib_name = ['libsim_model.so'];
				system(['gcc -O2 -fPIC -shared ', c_sources, ' -o ', lib_name]);
			end

			obj.C_sim_ext_fun = sim_create_ext_fun();

			% compile mex with model dependency
			if (strcmp(obj.opts_struct.compile_mex, 'true'))
				sim_compile_mex_model_dep(obj.model_struct, obj.opts_struct);
			end

			% get pointers for external functions in model
			if (strcmp(obj.opts_struct.method, 'erk'))
				obj.C_sim_ext_fun = sim_set_ext_fun_dyn_expl(obj.C_sim_ext_fun, obj.model_struct, obj.opts_struct);
			elseif (strcmp(obj.opts_struct.method, 'irk'))
				obj.C_sim_ext_fun = sim_set_ext_fun_dyn_impl(obj.C_sim_ext_fun, obj.model_struct, obj.opts_struct);
			else
				fprintf('\ncodegen_model: method not supported: %s\n', obj.opts_struct.method);
				return;
			end
			% set in model ( = casadi functions )
			sim_set_model(obj.opts_struct, obj.C_sim, obj.C_sim_ext_fun);
		end


		function set(obj, field, value)
			sim_set(obj.model_struct, obj.opts_struct, obj.C_sim, obj.C_sim_ext_fun, field, value);
		end


		function solve(obj)
			sim_solve(obj.C_sim);
		end


		function value = get(obj, field)
			value = sim_get(obj.C_sim, field);
		end


		function delete(obj)
			sim_destroy(obj.C_sim);
			sim_destroy_ext_fun(obj.model_struct, obj.C_sim_ext_fun);
		end


	end % methods



end % class

