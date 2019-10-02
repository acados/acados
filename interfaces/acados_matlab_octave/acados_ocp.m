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

classdef acados_ocp < handle

    properties
        C_ocp
        C_ocp_ext_fun
        model_struct
        opts_struct
		acados_ocp_nlp_json
    end % properties



    methods


        function obj = acados_ocp(model, opts)
            obj.model_struct = model.model_struct;
            obj.opts_struct = opts.opts_struct;

			% TODO(andrea): this is temporary. later on the solver_config
			% object will separate from the OCP object
			
			model.acados_ocp_nlp_json.solver_config.qp_solver = upper(obj.opts_struct.qp_solver);
			model.acados_ocp_nlp_json.solver_config.integrator_type = upper(obj.opts_struct.sim_method);
			model.acados_ocp_nlp_json.solver_config.nlp_solver_type = upper(obj.opts_struct.nlp_solver);
			model.acados_ocp_nlp_json.dims.N = upper(obj.opts_struct.param_scheme_N);
			obj.acados_ocp_nlp_json = model.acados_ocp_nlp_json;

            [~,~] = mkdir(obj.opts_struct.output_dir);
            addpath(obj.opts_struct.output_dir);

            % detect GNSF structure
            if (strcmp(obj.opts_struct.sim_method, 'irk_gnsf'))
                if (strcmp(obj.opts_struct.gnsf_detect_struct, 'true'))
                    obj.model_struct = detect_gnsf_structure(obj.model_struct);
                    generate_get_gnsf_structure(obj.model_struct, obj.opts_struct);
                else
                    obj.model_struct = get_gnsf_structure(obj.model_struct);
                end
            end

            % check if mex interface exists already
            if is_octave()
                mex_exists = exist( fullfile(obj.opts_struct.output_dir,...
                    '/ocp_create.mex'), 'file');
            else
                mex_exists = exist( fullfile(obj.opts_struct.output_dir,...
                    '/ocp_create.mexa64'), 'file');
            end

            % compile mex interface (without model dependency)
            if (strcmp(obj.opts_struct.compile_interface, 'true') || ~mex_exists)
                ocp_compile_interface(obj.opts_struct);
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



		function generate_c_code(obj)
			% generate C code for CasADi functions
			if (strcmp(obj.model_struct.dyn_type, 'explicit'))
				acados_template_mex.generate_c_code_explicit_ode(obj.acados_ocp_nlp_json.model);
			elseif (strcmp(obj.model_struct.dyn_type, 'implicit'))
				if (strcmp(obj.opts_struct.sim_method, 'irk'))
					opts.generate_hess = 1;
					acados_template_mex.generate_c_code_implicit_ode(obj.acados_ocp_nlp_json.model, opts);
				end
			end

			
			% set include and lib path
			acados_folder = getenv('ACADOS_INSTALL_DIR');
			obj.acados_ocp_nlp_json.acados_include_path = [acados_folder, '/include'];
			obj.acados_ocp_nlp_json.acados_lib_path = [acados_folder, '/lib'];
			% strip non-numerical data
			
			% model
			model.name = obj.acados_ocp_nlp_json.model.name;
			
			obj.acados_ocp_nlp_json.model = [];
			obj.acados_ocp_nlp_json.model = model;
			
			con_h.name = obj.acados_ocp_nlp_json.con_h.name;
			
			obj.acados_ocp_nlp_json.con_h = [];
			obj.acados_ocp_nlp_json.con_h = con_h;
			
			con_h_e.name = obj.acados_ocp_nlp_json.con_h_e.name;
			
			obj.acados_ocp_nlp_json.con_h_e = [];
			obj.acados_ocp_nlp_json.con_h_e = con_h_e;
			
			con_p.name = obj.acados_ocp_nlp_json.con_p.name;
			
			obj.acados_ocp_nlp_json.con_p = [];
			obj.acados_ocp_nlp_json.con_p = con_p;
			
			con_p_e.name = obj.acados_ocp_nlp_json.con_p_e.name;
			
			obj.acados_ocp_nlp_json.con_p_e = [];
			obj.acados_ocp_nlp_json.con_p_e = con_p_e;
			
			% post process numerical data (mostly cast scalars to 1-dimensional cells)
			constr = obj.acados_ocp_nlp_json.constraints;
			%props = properties(constr);
			props = fieldnames(constr);
			for iprop = 1:length(props)
				thisprop = props{iprop};
				%%%Add logic here if you want to work with select properties
				thisprop_value = constr.(thisprop);
				%%%Add logic here if you want to do something based on the property's value
				if size(thisprop_value) == [1 1]
					constr.(thisprop) = num2cell(constr.(thisprop));
				end
			end
			obj.acados_ocp_nlp_json.constraints = constr;
			
			cost = obj.acados_ocp_nlp_json.cost;
			%props = properties(cost);
			props = fieldnames(cost);
			for iprop = 1:length(props)
				thisprop = props{iprop};
				%%%Add logic here if you want to work with select properties
				thisprop_value = cost.(thisprop);
				%%%Add logic here if you want to do something based on the property's value
				if norm(size(thisprop_value) - [1, 1]) == 0
					cost.(thisprop) = num2cell(cost.(thisprop));
				end
			end
			obj.acados_ocp_nlp_json.cost = cost;
			
			% load JSON layout
			acados_folder = getenv('ACADOS_INSTALL_DIR');

			acados_layout = jsondecode(fileread([acados_folder,...
                '/interfaces/acados_template/acados_template/acados_layout.json']));

			dims = obj.acados_ocp_nlp_json.dims;
			% reshape constraints
			constr = obj.acados_ocp_nlp_json.constraints;
			constr_l = acados_layout.constraints;
			fields = fieldnames(constr_l);
			for i = 1:numel(fields)
				if strcmp(constr_l.(fields{i}){1}, 'ndarray')
					if length(constr_l.(fields{i}){2}) == 1
						this_dims = [dims.(constr_l.(fields{i}){2}{1}), 1];
					else
						this_dims = [dims.(constr_l.(fields{i}){2}{1}), dims.(constr_l.(fields{i}){2}{1})];
					end
					constr.(fields{i}) = reshape(constr.(fields{i}), this_dims);
				end
			end
			obj.acados_ocp_nlp_json.constraints = constr;
			
			% reshape cost
			cost = obj.acados_ocp_nlp_json.cost;
			cost_l = acados_layout.cost;
			fields = fieldnames(cost_l);
			for i = 1:numel(fields)
				if strcmp(cost_l.(fields{i}){1}, 'ndarray')
					if length(cost_l.(fields{i}){2}) == 1
						this_dims = [dims.(cost_l.(fields{i}){2}{1}), 1];
					else
						this_dims = [dims.(cost_l.(fields{i}){2}{1}), dims.(cost_l.(fields{i}){2}{2})];
					end
					cost.(fields{i}) = reshape(cost.(fields{i}), this_dims);
					% convert 1-dimensional arrays to cells
					if length(cost_l.(fields{i}){2}) == 2 && (this_dims(1) == 1 || this_dims(2) == 1)
						field_as_cell = {};
						for j = 1:max(this_dims(1), this_dims(2))
							field_as_cell{end+1} = num2cell(cost.(fields{i})(j));
						end
						cost.(fields{i}) = field_as_cell;
					end
				end
			end
			obj.acados_ocp_nlp_json.cost = cost;
			
			% dump JSON file
			json_string = jsonencode(obj.acados_ocp_nlp_json);
			fid = fopen('acados_ocp_nlp.json', 'w');
			if fid == -1, error('Cannot create JSON file'); end
			fwrite(fid, json_string, 'char');
			fclose(fid);
			% render templated C code
			% old call (Python + Jinja)
            % acados_template_mex.generate_solver('acados_ocp_nlp.json', '/home/andrea/.acados_t/bin/python3')
            acados_template_mex.generate_solver_matlab('acados_ocp_nlp.json')
		end




        function eval_param_sens(obj, field, stage, index)
            ocp_eval_param_sens(obj.C_ocp, field, stage, index);
        end


        function set(varargin)
            if nargin==3
                obj = varargin{1};
                field = varargin{2};
                value = varargin{3};
                ocp_set(obj.model_struct, obj.opts_struct, obj.C_ocp, obj.C_ocp_ext_fun, field, value);
            elseif nargin==4
                obj = varargin{1};
                field = varargin{2};
                value = varargin{3};
                stage = varargin{4};
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



        function print(varargin)
            if nargin < 2
                field = 'stat';
            else
                field = varargin{2};
            end

            obj = varargin{1};
            ocp_solver_string = obj.opts_struct.nlp_solver;

            if strcmp(field, 'stat')
                stat = obj.get('stat');
                if strcmp(ocp_solver_string, 'sqp')
                    fprintf('\niter\tres_stat\tres_eq\t\tres_ineq\tres_comp\tqp_stat\tqp_iter');
                    if size(stat,2)>7
                        fprintf('\tqp_res_stat\tqp_res_eq\tqp_res_ineq\tqp_res_comp');
                    end
                    fprintf('\n');
                    for jj=1:size(stat,1)
                        fprintf('%d\t%e\t%e\t%e\t%e\t%d\t%d', stat(jj,1), stat(jj,2), stat(jj,3), stat(jj,4), stat(jj,5), stat(jj,6), stat(jj,7));
                        if size(stat,2)>7
                            fprintf('\t%e\t%e\t%e\t%e', stat(jj,8), stat(jj,9), stat(jj,10), stat(jj,11));
                        end
                        fprintf('\n');
                    end
                    fprintf('\n');
                elseif strcmp(ocp_solver_string, 'sqp_rti')
                    fprintf('\niter\tqp_status\tqp_iter\n');
                    for jj=1:size(stat,1)
                        fprintf('%d\t%d\t\t%d', stat(jj,1), stat(jj,2), stat(jj,3));
                        fprintf('\n');
                    end
                end

            else
                fprintf('unsupported field in function print of acados_ocp, got %s', field);
                keyboard
            end

        end


        function delete(obj)
            if ~isempty(obj.C_ocp_ext_fun)
                ocp_destroy_ext_fun(obj.model_struct, obj.C_ocp, obj.C_ocp_ext_fun);
            end
            if ~isempty(obj.C_ocp) 
                ocp_destroy(obj.C_ocp);
            end
        end


    end % methods



end % class

