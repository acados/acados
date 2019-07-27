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

			[~,~] = mkdir(obj.opts_struct.output_dir);
			addpath(obj.opts_struct.output_dir);

			% detect GNSF structure
			if (strcmp(obj.opts_struct.method, 'irk_gnsf'))
				if (strcmp(obj.opts_struct.gnsf_detect_struct, 'true'))
					obj.model_struct = detect_gnsf_structure(obj.model_struct);
					generate_get_gnsf_structure(obj.model_struct, obj.opts_struct);
				else
					obj.model_struct = get_gnsf_structure(obj.model_struct);
				end
			end

			% compile mex without model dependency
			if (strcmp(obj.opts_struct.compile_mex, 'true'))
				sim_compile_mex(obj.opts_struct);
			end

			obj.C_sim = sim_create(obj.model_struct, obj.opts_struct);

			% generate and compile casadi functions
			if (strcmp(obj.opts_struct.codgen_model, 'true'))
				sim_generate_casadi_ext_fun(obj.model_struct, obj.opts_struct)
			end

			obj.C_sim_ext_fun = sim_create_ext_fun();

			% compile mex with model dependency & set pointers for external functions in model
			obj.C_sim_ext_fun = sim_set_ext_fun(obj.C_sim, obj.C_sim_ext_fun, obj.model_struct, obj.opts_struct);

			% precompute
			sim_precompute(obj.C_sim);

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
			if ~isempty(obj.C_sim)
				sim_destroy(obj.C_sim);
			end
			if ~isempty(obj.C_sim_ext_fun)
				sim_destroy_ext_fun(obj.model_struct, obj.C_sim_ext_fun);
			end
		end


	end % methods



end % class

