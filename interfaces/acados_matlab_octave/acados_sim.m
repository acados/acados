%
% Copyright (c) The acados authors.
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
        % templated solver
        t_sim
        % matlab objects
        model_struct
        opts_struct
        acados_sim_json
        code_gen_dir
    end % properties



    methods


        function obj = acados_sim(model, opts)
            obj.model_struct = model.model_struct;
            obj.opts_struct = opts.opts_struct;

            [~,~] = mkdir(obj.opts_struct.output_dir);
            addpath(obj.opts_struct.output_dir);

            % check model consistency
            obj.model_struct = create_consistent_empty_fields(obj.model_struct, obj.opts_struct);
            % detect GNSF structure
            if (strcmp(obj.opts_struct.method, 'irk_gnsf'))
                if (strcmp(obj.opts_struct.gnsf_detect_struct, 'true'))
                    obj.model_struct = detect_gnsf_structure(obj.model_struct);
                    generate_get_gnsf_structure(obj.model_struct, obj.opts_struct);
                else
                    obj.model_struct = get_gnsf_structure(obj.model_struct);
                end
            end

            % check if path contains spaces
            if ~isempty(strfind(obj.opts_struct.output_dir, ' '))
                error(strcat('acados_ocp: Path should not contain spaces, got: ',...
                    obj.opts_struct.output_dir));
            end

            %% compile mex without model dependency
            % check if mex interface exists already
            if strcmp(obj.opts_struct.compile_interface, 'true')
                compile_interface = true;
            elseif strcmp(obj.opts_struct.compile_interface, 'false')
                compile_interface = false;
            elseif strcmp(obj.opts_struct.compile_interface, 'auto')
                if is_octave()
                    extension = '.mex';
                else
                    extension = ['.' mexext];
                end
                compile_interface = ~exist(fullfile(obj.opts_struct.output_dir, ['/sim_create', extension]), 'file');
            else
                obj.model_struct.cost_type
                error('acados_sim: field compile_interface is , supported values are: true, false, auto');
            end

            if ( compile_interface )
                sim_compile_interface(obj.opts_struct);
            end

            % detect dimensions & sanity checks
            obj.model_struct = detect_dims_sim(obj.model_struct,obj.opts_struct);


            % create template sim
            obj.acados_sim_json = set_up_acados_sim_json(obj);
            sim_generate_c_code(obj);

            % templated MEX
            return_dir = pwd();
            obj.code_gen_dir = obj.acados_sim_json.code_export_directory; 
            cd(obj.code_gen_dir)

            mex_sim_solver = str2func(sprintf('%s_mex_sim_solver', obj.model_struct.name));
            obj.t_sim = mex_sim_solver();
            addpath(pwd());

            cd(return_dir)
        end


        function set(obj, field, value)
            obj.t_sim.set(field, value);
        end


        function status = solve(obj)
            status = obj.t_sim.solve();
        end


        function value = get(obj, field)
            value = obj.t_sim.get(field);
        end


        % function delete(obj)
        %     Use default implementation.
        %     MATLAB destroys the property values after the destruction of the object. 
        %     Because `t_sim` is the only referrence to the `mex_sim_solver` object, MATLAB also destroys the latter.
        % end


    end % methods



end % class

