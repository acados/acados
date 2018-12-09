classdef acados_integrator < handle
	


	properties
		C_sim
		model
		opts
	end % properties



	methods


		function obj = acados_integrator(model, opts)
			obj.model = model;
			obj.opts = opts;
			% create opts struct
			opts_struct = struct;
			opts_struct.num_stages = opts.num_stages;
			opts_struct.num_steps = opts.num_steps;
			if (strcmp(opts.sens_forw, 'true'))
				opts_struct.sens_forw = 1;
			else
				opts_struct.sens_forw = 0;
			end
			opts_struct.scheme = opts.scheme;
			obj.C_sim = sim_create(model, opts_struct);
		end


		function codegen_model(obj)
			if obj.opts.codgen_model
				if obj.opts.scheme=='erk'
					% generate c for function and derivatives using casadi
					generate_c_code_explicit_ode(obj.model);
					% compile the code in a shared library
					c_sources = ' ';
					c_sources = [c_sources, obj.model.name, '_expl_ode_fun.c ']; 
					c_sources = [c_sources, obj.model.name, '_expl_vde_forw.c ']; 
					c_sources = [c_sources, obj.model.name, '_expl_vde_adj.c ']; 
					c_sources = [c_sources, obj.model.name, '_expl_ode_hess.c ']; 
					lib_name = [obj.model.name, '_expl.so'];
					system(['gcc -fPIC -shared ', c_sources, ' -o ', lib_name]);
				else
					fprintf('\nscheme not supported: %s\n', obj.opts.scheme);
				end
			end
		end


		function set(obj, field, value)
		end

	
		function delete(obj)
%			fprintf('\nin delete\n');
			sim_destroy(obj.C_sim)
		end


	end % methods



end % class





