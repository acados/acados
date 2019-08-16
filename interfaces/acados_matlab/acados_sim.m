%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;
%

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

