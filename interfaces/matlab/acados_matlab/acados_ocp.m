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


		function delete(obj)
%			fprintf('\nin delete\n');
			ocp_destroy(obj.C_ocp);
%			ocp_ext_fun_destroy(obj.opts_struct, obj.C_sim_ext_fun);
		end


	end % methods



end % class

