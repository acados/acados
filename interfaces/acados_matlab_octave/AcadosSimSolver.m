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

    properties (Access = public)
        sim % MATLAB class AcadosSim describing the initial value problem
    end

    properties (Access = private)
        t_sim % templated solver
        name
    end % properties

    methods

        function obj = AcadosSimSolver(sim, varargin)
            %% optional arguments:
            % varargin{1}: solver_creation_opts: this is a struct in which some of the fields can be defined to overwrite the default values.
            % The fields are:
            % - json_file: path to the json file containing the ocp description
            % - build: boolean, if true, the problem specific shared library is compiled
            % - generate: boolean, if true, the C code is generated
            % - compile_mex_wrapper: boolean, if true, the mex wrapper is compiled
            % - compile_interface: can be [], true or false. If [], the interface is compiled if it does not exist.
            % - output_dir: path to the directory where the MEX interface is compiled
            obj.sim = sim;

            % optional arguments
            % solver creation options
            default_solver_creation_opts = struct('json_file', '', ...
                    'build', true, ...
                    'generate', true, ...
                    'compile_mex_wrapper', true, ...
                    'compile_interface', [], ...
                    'output_dir', fullfile(pwd, 'build'));
            if length(varargin) > 0
                solver_creation_opts = varargin{1};
                % set non-specified opts to default
                fields = fieldnames(default_solver_creation_opts);
                for i = 1:length(fields)
                    if ~isfield(solver_creation_opts, fields{i})
                        solver_creation_opts.(fields{i}) = default_solver_creation_opts.(fields{i});
                    end
                end
            else
                solver_creation_opts = default_solver_creation_opts;
            end

            if isempty(sim) && isempty(solver_creation_opts.json_file)
                error('AcadosSimSolver: provide either a sim object or a json file');
            end

            if isempty(sim)
                json_file = solver_creation_opts.json_file;
            else
                % formulation provided
                if ~isempty(solver_creation_opts.json_file)
                    sim.json_file = solver_creation_opts.json_file;
                end
                json_file = sim.json_file;
                if ~isempty(sim.solver_options.compile_interface) && ~isempty(solver_creation_opts.compile_interface)
                    error('AcadosOcpSolver: provide either compile_interface in OCP object or solver_creation_opts');
                end
                if ~isempty(sim.solver_options.compile_interface)
                    solver_creation_opts.compile_interface = sim.solver_options.compile_interface;
                end
                % make consistent
                sim.make_consistent();
            end

            % compile mex sim interface if needed
            obj.compile_mex_sim_interface_if_needed(solver_creation_opts);

            %% generate
            if solver_creation_opts.generate
                obj.generate();
            end

            % load json: TODO!?
            acados_folder = getenv('ACADOS_INSTALL_DIR');
            addpath(fullfile(acados_folder, 'external', 'jsonlab'));
            acados_sim_struct = loadjson(fileread(json_file), 'SimplifyCell', 0);
            obj.name = acados_sim_struct.model.name;
            code_export_directory = acados_sim_struct.code_export_directory;

            %% compile problem specific shared library
            if solver_creation_opts.build
                obj.compile_sim_shared_lib(code_export_directory);
            end

            %% create solver
            return_dir = pwd();
            cd(code_export_directory)

            mex_sim_solver = str2func(sprintf('%s_mex_sim_solver', obj.name));
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


        function x_next = simulate(obj, x, u, z, xdot, p)
        % Simulate the system forward for the given x, u, p and return x_next.
        % The values xdot, z are used as initial guesses for implicit integrators, if provided.
        % Wrapper around solve() taking care of setting/getting inputs/outputs.
        % Fields which are set to an empty array or not provided will not be set.

            if nargin >= 2 && ~isempty(x)
                obj.set('x', x);
            end
            if nargin >= 3 && ~isempty(u)
                obj.set('u', u);
            end

            if strcmp(obj.sim.solver_options.integrator_type, 'IRK')
                if nargin >= 4 && ~isempty(z)
                    obj.set('z', z);
                end
                if nargin >= 5 && ~isempty(xdot)
                    obj.set('xdot', xdot);
                end
            end

            if nargin >= 6 && ~isempty(p)
                obj.set('p', p);
            end

            status = obj.solve();

            if status ~= 0
                error('AcadosSimSolver for model %s returned status %d.', obj.name, status);
            end

            x_next = obj.get('xn');
        end


        % function delete(obj)
        %     Use default implementation.
        %     MATLAB destroys the property values after the destruction of the object.
        %     Because `t_sim` is the only referrence to the `mex_sim_solver` object, MATLAB also destroys the latter.
        % end

    end % methods

    methods (Access = private)
        function generate(obj)
            % generate
            check_dir_and_create(fullfile(pwd, obj.sim.code_export_directory));
            obj.sim.generate_external_functions();

            obj.sim.dump_to_json()
            obj.sim.render_templates()
        end

        function compile_mex_sim_interface_if_needed(obj, solver_creation_opts)

            [~,~] = mkdir(solver_creation_opts.output_dir);
            addpath(solver_creation_opts.output_dir);

            % check if path contains spaces
            if ~isempty(strfind(solver_creation_opts.output_dir, ' '))
                error(strcat('compile_mex_sim_interface_if_needed: Path should not contain spaces, got: ',...
                    solver_creation_opts.output_dir));
            end

            %% compile mex without model dependency
            % check if mex interface exists already
            if isempty(solver_creation_opts.compile_interface) % auto-detect
                if is_octave()
                    extension = '.mex';
                else
                    extension = ['.' mexext];
                end
                solver_creation_opts.compile_interface = ~exist(fullfile(solver_creation_opts.output_dir, ['/sim_create', extension]), 'file');
            end

            if solver_creation_opts.compile_interface
                sim_compile_interface(solver_creation_opts.output_dir);
            end
        end

        function compile_sim_shared_lib(obj, export_dir)
            return_dir = pwd;
            cd(export_dir);
            if isunix
                [ status, result ] = system('make sim_shared_lib');
                if status
                    cd(return_dir);
                    error('Building templated code as shared library failed.\nGot status %d, result: %s',...
                        status, result);
                end
            else
                % check compiler
                use_msvc = false;
                if ~is_octave()
                    mexOpts = mex.getCompilerConfigurations('C', 'Selected');
                    if contains(mexOpts.ShortName, 'MSVC')
                        use_msvc = true;
                    end
                end
                % compile on Windows platform
                if use_msvc
                    % get env vars for MSVC
                    % msvc_env = fullfile(mexOpts.Location, 'VC\Auxiliary\Build\vcvars64.bat');
                    % assert(isfile(msvc_env), 'Cannot find definition of MSVC env vars.');
                    % detect MSVC version
                    msvc_ver_str = "Visual Studio " + mexOpts.Version(1:2) + " " + mexOpts.Name(22:25);
                    [ status, result ] = system(['cmake -G "' + msvc_ver_str + '" -A x64 -DCMAKE_BUILD_TYPE=Release -DBUILD_ACADOS_SIM_SOLVER_LIB=ON -DBUILD_ACADOS_OCP_SOLVER_LIB=OFF -S . -B .']);
                else
                    [ status, result ] = system('cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DBUILD_ACADOS_SIM_SOLVER_LIB=ON -DBUILD_ACADOS_OCP_SOLVER_LIB=OFF -S . -B .');
                end
                if status
                    cd(return_dir);
                    error('Generating buildsystem failed.\nGot status %d, result: %s',...
                        status, result);
                end
                [ status, result ] = system('cmake --build . --config Release');
                if status
                    cd(return_dir);
                    error('Building templated code as shared library failed.\nGot status %d, result: %s',...
                        status, result);
                end
            end

            cd(return_dir);
        end % methods (Access = private)

    end
end % class

