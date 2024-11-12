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
        t_sim % templated solver
        sim % Matlab class AcadosSim describing the initial value problem
    end % properties



    methods

        function obj = AcadosSimSolver(sim, output_dir)

            if nargin < 2
                output_dir = fullfile(pwd, 'build');
            end

            % check model consistency
            obj.sim = sim;
            sim.make_consistent();

            % create template sim
            sim_generate_c_code(obj.sim);

            % compile mex sim interface if needed
            obj.compile_mex_sim_interface_if_needed(output_dir)

            % templated MEX
            return_dir = pwd();
            cd(obj.sim.code_export_directory)

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

    methods (Access = private)
        function compile_mex_sim_interface_if_needed(obj, output_dir)

            [~,~] = mkdir(output_dir);
            addpath(output_dir);

            % check if path contains spaces
            if ~isempty(strfind(output_dir, ' '))
                error(strcat('acados_ocp: Path should not contain spaces, got: ',...
                    output_dir));
            end

            %% compile mex without model dependency
            % check if mex interface exists already
            if isempty(obj.sim.solver_options.compile_interface) % auto-detect
                if is_octave()
                    extension = '.mex';
                else
                    extension = ['.' mexext];
                end
                obj.sim.solver_options.compile_interface = ~exist(fullfile(output_dir, ['/sim_create', extension]), 'file');
            end

            if obj.sim.solver_options.compile_interface
                sim_compile_interface(output_dir);
            end
        end
    end
end % class

