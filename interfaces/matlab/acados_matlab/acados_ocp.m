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
			obj.C_ocp = ocp_create(obj.model_struct, obj.opts_struct);
		end



		function codegen_model(obj)
			if obj.opts_struct.codgen_model
				c_sources = ' ';
				if (strcmp(obj.opts_struct.sim_solver, 'erk'))
					% generate c for function and derivatives using casadi
					generate_c_code_explicit_ode(obj.model_struct);
					% compile the code in a shared library
					c_sources = [c_sources, 'model_expl_ode_fun.c '];
					c_sources = [c_sources, 'model_expl_vde_for.c '];
					c_sources = [c_sources, 'model_expl_vde_adj.c '];
					c_sources = [c_sources, 'model_expl_ode_hes.c '];
				elseif (strcmp(obj.opts_struct.sim_solver, 'irk'))
					% generate c for function and derivatives using casadi
					generate_c_code_implicit_ode(obj.model_struct);
					% compile the code in a shared library
					c_sources = [c_sources, 'model_impl_ode_fun.c '];
					c_sources = [c_sources, 'model_impl_ode_fun_jac_x_xdot.c '];
					c_sources = [c_sources, 'model_impl_ode_fun_jac_x_xdot_u.c '];
					c_sources = [c_sources, 'model_impl_ode_jac_x_xdot_u.c '];
				else
					fprintf('\ncodegen_model: sim solver not supported: %s\n', obj.opts_struct.sim_solver);
					return;
				end
				lib_name = ['libmodel.so'];
				system(['gcc -fPIC -shared ', c_sources, ' -o ', lib_name]);
			end

			% get acados folder (if set)
			acados_folder = getenv('ACADOS_FOLDER');
			% default folder
			if length(acados_folder) == 0
				acados_folder = '../../../';
			end
			% set paths
			acados_include = ['-I' acados_folder];
			acados_interfaces_include = ['-I' acados_folder, 'interfaces'];
			acados_lib_path = ['-L' acados_folder, 'lib'];
			acados_matlab_lib_path = ['-L' acados_folder, 'interfaces/matlab/acados_matlab/'];

			% get pointers for external functions in model
			if (strcmp(obj.opts_struct.sim_solver, 'erk'))
				mex(acados_include, acados_interfaces_include, acados_lib_path, acados_matlab_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', '-lmodel', 'ocp_expl_ext_fun_create.c');
				obj.C_ocp_ext_fun = ocp_expl_ext_fun_create(obj.opts_struct);
			elseif (strcmp(obj.opts_struct.sim_solver, 'irk'))
				mex(acados_include, acados_interfaces_include, acados_lib_path, acados_matlab_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', '-lmodel', 'ocp_impl_ext_fun_create.c');
				obj.C_ocp_ext_fun = ocp_impl_ext_fun_create(obj.opts_struct);
			else
				fprintf('\ncodegen_model: sim_solver not supported: %s\n', obj.opts_struct.sim_solver);
				return;
			end
			% set in model ( = casadi functions )
			ocp_set_model(obj.opts_struct, obj.C_ocp, obj.C_ocp_ext_fun);
		end



		function solve(obj)
			ocp_solve(obj.C_ocp);
		end



		function get(obj, field, value)
			ocp_get(obj.C_ocp, field, value);
		end



		function delete(obj)
%			fprintf('\nin delete\n');
			ocp_destroy(obj.C_ocp);
			ocp_ext_fun_destroy(obj.opts_struct, obj.C_ocp_ext_fun);
		end


	end % methods



end % class

