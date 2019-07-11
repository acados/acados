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

			% detect GNSF structure
			if (strcmp(obj.opts_struct.method, 'irk_gnsf'))
				obj.model_struct = detect_gnsf_structure(obj.model_struct);
				obj.model_struct
				% TODO code-generate a file with e.g. dims, ..., to avoid to detect the gnsf structure all the time
			end

			% compile mex without model dependency
			if (strcmp(obj.opts_struct.compile_mex, 'true'))
				sim_compile_mex();
			end

			obj.C_sim = sim_create(obj.model_struct, obj.opts_struct);

			% generate and compile casadi functions
			if (strcmp(obj.opts_struct.codgen_model, 'true'))
				sim_compile_casadi_functions(obj.model_struct, obj.opts_struct)
			end

			obj.C_sim_ext_fun = sim_create_ext_fun();

			% compile mex with model dependency
			if (strcmp(obj.opts_struct.compile_mex, 'true'))
				sim_compile_mex_model_dep(obj.model_struct, obj.opts_struct);
			end

			% get pointers for external functions in model
			if (strcmp(obj.opts_struct.method, 'erk'))
				obj.C_sim_ext_fun = sim_set_ext_fun_dyn_expl(obj.C_sim, obj.C_sim_ext_fun, obj.model_struct, obj.opts_struct);
			elseif (strcmp(obj.opts_struct.method, 'irk'))
				obj.C_sim_ext_fun = sim_set_ext_fun_dyn_impl(obj.C_sim, obj.C_sim_ext_fun, obj.model_struct, obj.opts_struct);
			else
				fprintf('\ncodegen_model: method not supported: %s\n', obj.opts_struct.method);
				return;
			end
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

