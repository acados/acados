classdef acados_integrator < handle

	properties
		C_sim
		C_sim_ext_fun
		model
		opts
		opts_struct
	end % properties



	methods


		function obj = acados_integrator(model, opts)
			obj.model = model;
			obj.opts = opts;
			% create opts struct
			obj.opts_struct = struct;
			obj.opts_struct.num_stages = opts.num_stages;
			obj.opts_struct.num_steps = opts.num_steps;
			if (strcmp(opts.sens_forw, 'true'))
				obj.opts_struct.sens_forw = 1;
			else
				obj.opts_struct.sens_forw = 0;
			end
			obj.opts_struct.scheme = opts.scheme;
			obj.C_sim = sim_create(model, obj.opts_struct);
		end


		function codegen_model(obj)
			if obj.opts.codgen_model
				c_sources = ' ';
				if (strcmp(obj.opts.scheme, 'erk'))
					% generate c for function and derivatives using casadi
					generate_c_code_explicit_ode(obj.model);
					% compile the code in a shared library
					c_sources = [c_sources, 'model_expl_ode_fun.c '];
					c_sources = [c_sources, 'model_expl_vde_for.c '];
					c_sources = [c_sources, 'model_expl_vde_adj.c '];
					c_sources = [c_sources, 'model_expl_ode_hes.c '];
				elseif (strcmp(obj.opts.scheme, 'irk'))
					% generate c for function and derivatives using casadi
					generate_c_code_implicit_ode(obj.model);
					% compile the code in a shared library
					c_sources = [c_sources, 'model_impl_ode_fun.c '];
					c_sources = [c_sources, 'model_impl_ode_fun_jac_x_xdot.c '];
					c_sources = [c_sources, 'model_impl_ode_fun_jac_x_xdot_u.c '];
					c_sources = [c_sources, 'model_impl_ode_jac_x_xdot_u.c '];
				else
					fprintf('\ncodegen_model: scheme not supported: %s\n', obj.opts.scheme);
					return;
				end
				lib_name = ['libmodel.so'];
				system(['gcc -fPIC -shared ', c_sources, ' -o ', lib_name]);

				acados_folder = getenv('ACADOS_FOLDER');

				if length(acados_folder) == 0
					acados_folder = '../../../';
				end

				include_acados = ['-I' acados_folder];
				include_interfaces = ['-I' acados_folder, 'interfaces'];
				acados_lib_path = ['-L' acados_folder, 'lib'];
				acados_matlab_lib_path = ['-L' acados_folder, 'interfaces/matlab/acados_matlab/'];


				% get pointers for external functions in model
				if (strcmp(obj.opts.scheme, 'erk'))
					mex(include_acados, include_interfaces, acados_lib_path, acados_matlab_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', '-lmodel', 'sim_expl_ext_fun_create.c');
					obj.C_sim_ext_fun = sim_expl_ext_fun_create(obj.opts_struct);
				elseif (strcmp(obj.opts.scheme, 'irk'))
					mex(include_acados, include_interfaces, acados_lib_path, acados_matlab_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', '-lmodel', 'sim_impl_ext_fun_create.c');
					obj.C_sim_ext_fun = sim_impl_ext_fun_create(obj.opts_struct);
				else
					fprintf('\ncodegen_model: scheme not supported: %s\n', obj.opts.scheme);
					return;
				end
				% set in model ( = casadi functions )
				sim_set_model(obj.opts_struct, obj.C_sim, obj.C_sim_ext_fun);
			end
		end


		function set(obj, field, value)
			sim_set(obj.C_sim, field, value);
		end


		function solve(obj)
			sim_solve(obj.C_sim);
		end


		function get(obj, field, value)
			sim_get(obj.C_sim, field, value);
		end

		function delete(obj)
%			fprintf('\nin delete\n');
			sim_destroy(obj.C_sim);
			sim_ext_fun_destroy(obj.opts_struct, obj.C_sim_ext_fun);
		end


	end % methods



end % class





