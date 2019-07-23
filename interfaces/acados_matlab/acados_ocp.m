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

			% create build folder and add to path
      addpath('build');
      rmpath('build')
			[~] = rmdir('build', 's');
      [~,~] = mkdir('build');
			addpath('build');

			% detect GNSF structure
			if (strcmp(obj.opts_struct.sim_method, 'irk_gnsf'))
				if (strcmp(obj.opts_struct.gnsf_detect_struct, 'true'))
					obj.model_struct = detect_gnsf_structure(obj.model_struct);
					generate_get_gnsf_structure(obj.model_struct);
				else
					obj.model_struct = get_gnsf_structure(obj.model_struct);
				end
			end

			% compile mex without model dependency
			if (strcmp(obj.opts_struct.compile_mex, 'true'))
				ocp_compile_mex(obj.opts_struct);
			end

			obj.C_ocp = ocp_create(obj.model_struct, obj.opts_struct);

			% generate and compile casadi functions
			if (strcmp(obj.opts_struct.codgen_model, 'true'))
				ocp_generate_casadi_ext_fun(obj.model_struct, obj.opts_struct);
			end

			obj.C_ocp_ext_fun = ocp_create_ext_fun();

			% compile mex with model dependency & set pointers for external functions in model
			obj.C_ocp_ext_fun = ocp_set_ext_fun(obj.C_ocp, obj.C_ocp_ext_fun, obj.model_struct, obj.opts_struct);

			% precompute
			ocp_precompute(obj.C_ocp);

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

