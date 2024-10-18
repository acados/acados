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

classdef AcadosSim < handle

    properties
        % file structure
        acados_include_path
        acados_lib_path
        shared_lib_ext
        json_file
        cython_include_dirs
        code_export_directory
        % struct / object
        dims
        model
        solver_options
        % plain data
        parameter_values
        problem_class
        external_function_files_model
    end

    methods
        function obj = AcadosSim()
            % most fields are initialized as a placeholder
            obj.acados_include_path = [];
            obj.acados_lib_path = [];
            obj.shared_lib_ext = [];
            obj.json_file = [];
            obj.cython_include_dirs = [];
            obj.code_export_directory = [];

            obj.dims = AcadosSimDims();
            obj.model = AcadosModel();

            obj.solver_options = AcadosSimOptions();

            obj.parameter_values = [];
            obj.problem_class = 'SIM';
        end

        function make_consistent(self)

            % file structures
            acados_folder = getenv('ACADOS_INSTALL_DIR');
            self.acados_include_path = [acados_folder, '/include'];
            self.acados_lib_path = [acados_folder, '/lib'];
            self.shared_lib_ext = '.so';
            if ismac()
                self.shared_lib_ext = '.dylib';
            end
            if isempty(self.json_file)
                self.json_file = 'acados_sim.json';
            end

            if isempty(self.code_export_directory)
                self.code_export_directory = fullfile(pwd(), 'c_generated_code');
            end

            % model
            self.model.make_consistent(self.dims);

            if self.dims.np_global > 0
                error('p_global is not supported for AcadosSim.')
            end

            % detect GNSF structure
            if strcmp(self.solver_options.integrator_type, 'GNSF')
                % TODO: interface these options
                gnsf_transcription_opts = struct();
                if self.dims.gnsf_nx1 + self.dims.gnsf_nx2 ~= self.dims.nx
                    detect_gnsf_structure(self.model, self.dims, gnsf_transcription_opts);
                else
                    warning('No GNSF model detected, assuming all required fields are set.')
                end
            end

            % parameters
            if self.dims.np > 0
                if isempty(self.parameter_values)
                    warning(['opts_struct.parameter_values are not set.', ...
                                10 'Using zeros(np,1) by default.' 10 'You can update them later using set().']);
                    self.parameter_values = zeros(self.dims.np,1);
                end
            end

            % solver options checks
            opts = self.solver_options;
            if ~strcmp(opts.integrator_type, 'ERK') && ~strcmp(opts.integrator_type, 'IRK') && ...
               ~strcmp(opts.integrator_type, 'GNSF') && ~strcmp(opts.integrator_type, 'DISCRETE')
                error(['integrator_type = ', opts.integrator_type, ' not available. Choose ERK, IRK, GNSF.']);
            end

            if length(opts.num_stages) ~= 1
                error('num_stages should be a scalar.');
            end
            if length(opts.num_steps) ~= 1
                error('num_steps should be a scalar.');
            end
            if length(opts.newton_iter) ~= 1
                error('newton_iter should be a scalar.');
            end
            % check bool options
            if ~islogical(opts.sens_forw)
                error('sens_forw should be a boolean.');
            end
            if ~islogical(opts.sens_adj)
                error('sens_adj should be a boolean.');
            end
            if ~islogical(opts.sens_algebraic)
                error('sens_algebraic should be a boolean.');
            end
            if ~islogical(opts.sens_hess)
                error('sens_hess should be a boolean.');
            end
            if ~islogical(opts.output_z)
                error('output_z should be a boolean.');
            end
            if ~strcmp(opts.collocation_type, "GAUSS_LEGENDRE") && ~strcmp(opts.collocation_type, "GAUSS_RADAU_IIA")
                error(['collocation_type = ', opts.collocation_type, ' not available. Choose GAUSS_LEGENDRE, GAUSS_RADAU_IIA.']);
            end

            if(strcmp(self.solver_options.integrator_type, "ERK"))
                if(self.solver_options.num_stages == 1 || self.solver_options.num_stages == 2 || ...
                    self.solver_options.num_stages == 3 || self.solver_options.num_stages == 4)
                else
                    error(['ERK: num_stages = ', num2str(self.solver_options.num_stages) ' not available. Only number of stages = {1,2,3,4} implemented!']);
                end
            end



        end

        function s = struct(self)
            if exist('properties')
                publicProperties = eval('properties(self)');
            else
                publicProperties = fieldnames(self);
            end
            s = struct();
            for fi = 1:numel(publicProperties)
                s.(publicProperties{fi}) = self.(publicProperties{fi});
            end
        end
    end

end % class
