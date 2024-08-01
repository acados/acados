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

classdef AcadosSimSolver < handle

    properties
        % templated solver
        t_sim
        % matlab objects
        sim
        code_gen_dir
    end % properties



    methods

        function obj = AcadosSimSolver(model, opts)

            % TODO where to get these from?
            output_dir = opts.opts_struct.output_dir;
            gnsf_transcription_opts = struct();

            [~,~] = mkdir(opts.opts_struct.output_dir);
            addpath(opts.opts_struct.output_dir);

            % check model consistency
            obj.sim = setup_sim(model.model_struct, opts.opts_struct);
            obj.sim.model.make_consistent(obj.sim.dims);
            % detect dimensions & sanity checks
            detect_dims_sim(obj.sim.model, obj.sim.sim_options);

            % detect GNSF structure
            if strcmp(obj.sim.sim_options.integrator_type, 'GNSF')
                if obj.sim.dims.gnsf_nx1 + obj.sim.dims.gnsf_nx2 ~= obj.sim.dims.nx
                    detect_gnsf_structure(obj.sim.model, obj.sim.dims, gnsf_transcription_opts);
                else
                    warning('No GNSF model detected, assuming all required fields are set.')
                end
            end

            % parameters
            if obj.sim.dims.np > 0
                if isempty(obj.sim.parameter_values)
                    warning(['opts_struct.parameter_values are not set.', ...
                                10 'Using zeros(np,1) by default.' 10 'You can update them later using set().']);
                    obj.sim.parameter_values = zeros(obj.sim.dims.np,1);
                end
            end

            % check if path contains spaces
            if ~isempty(strfind(output_dir, ' '))
                error(strcat('acados_ocp: Path should not contain spaces, got: ',...
                    output_dir));
            end

            %% compile mex without model dependency
            % check if mex interface exists already
            if strcmp(opts.opts_struct.compile_interface, 'true')
                compile_interface = true;
            elseif strcmp(opts.opts_struct.compile_interface, 'false')
                compile_interface = false;
            elseif strcmp(opts.opts_struct.compile_interface, 'auto')
                if is_octave()
                    extension = '.mex';
                else
                    extension = ['.' mexext];
                end
                compile_interface = ~exist(fullfile(output_dir, ['/sim_create', extension]), 'file');
            else
                error('acados_sim: field compile_interface is , supported values are: true, false, auto');
            end

            if (compile_interface)
                sim_compile_interface(output_dir);
            end

            % create template sim
            sim_generate_c_code(obj.sim);

            % templated MEX
            return_dir = pwd();
            obj.code_gen_dir = obj.sim.code_export_directory;
            cd(obj.code_gen_dir)

            mex_sim_solver = str2func(sprintf('%s_mex_sim_solver', obj.sim.model.name));
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

