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

classdef AcadosOcpSimulinkOptions < handle

    properties
        outputs
        inputs
        samplingtime
        show_port_info
        zoro_iterations
        generate_simulink_block
        customizable_inputs
    end

    methods
        function obj = AcadosOcpSimulinkOptions(problem_class)
            if nargin < 1
                problem_class = 'OCP';
            end
            if ~ismember(problem_class, {'OCP', 'MOCP'})
                error('problem_class must be ''OCP'' or ''MOCP''.');
            end
            obj.outputs = AcadosOcpSimulinkOutputs();
            obj.inputs = AcadosOcpSimulinkInputs();
            obj.samplingtime = 't0';
            obj.show_port_info = 1;
            obj.zoro_iterations = 1;
            obj.generate_simulink_block = 1;
            obj.customizable_inputs = struct();
            if strcmp(problem_class, 'MOCP')
                fields_to_disable = [AcadosOcpSimulinkOptions.nonsupported_mocp_inputs(), ...
                                    {'lbx', 'ubx', 'lbx_e', 'ubx_e', 'lh', 'uh', 'y_ref_0', 'y_ref_e'}];
                for i = 1:length(fields_to_disable)
                    obj.inputs.(fields_to_disable{i}) = 0;
                end
            end
        end

        function add_customizable_input(self, input_name, input_spec)
            self.customizable_inputs.(input_name) = input_spec;
        end

        function add_sparse_param_port(self, idx_p, port_name, stage_idx_0, stage_idx_e)
            %% add_sparse_param_port:
            % allows one to specify information for an input port of the simulink block corresponding to an acados OCP solver.
            %
            % inputs:
            % idx_p is a 0-based vector of parameter indices to be updated by the port.
            % port_name is used to identify the port and print information
            % stage_idx_0 is the first stage for which the parameters should be updated by the port
            % stage_idx_e is the last stage for which the parameters should be updated by the port

            % sanity checks:
            if stage_idx_0 > stage_idx_e
                error("stage_idx_0 > stage_idx_e")
            end

            input_name = strcat('sparse_parameter_', port_name);

            input_spec = struct('parameter_indices', idx_p, 'stage_idx_0', stage_idx_0, 'stage_idx_e', stage_idx_e);

            % NOTE: order matters here.
            % struct() auto-expands cell-array args into a struct array, so num2cell(idx_p)
            % can't be passed directly. Assign it after, and wrap scalars in a 1x1 cell so
            % jsonencode emits a JSON array instead of a bare number.
            if length(idx_p) == 1
                input_spec.parameter_indices = reshape(num2cell(idx_p), [1, 1]);
            end

            self.add_customizable_input(input_name, input_spec);
        end

        function s = to_struct(self)
            if exist('properties')
                publicProperties = eval('properties(self)');
            else
                publicProperties = fieldnames(self);
            end

            s = struct();
            for fi = 1:numel(publicProperties)
                f = publicProperties{fi};
                value = self.(f);
                if isa(value, 'AcadosOcpSimulinkInputs') || isa(value, 'AcadosOcpSimulinkOutputs')
                    s.(f) = value.to_struct();
                else
                    s.(f) = value;
                end
            end

            if isempty(fieldnames(self.customizable_inputs))
                s = rmfield(s, 'customizable_inputs');
            end
            s = orderfields(s);
        end

        function make_consistent(self, solver_options, problem_class)
            if self.inputs.rti_phase && solver_options.nlp_solver_type ~= 'SQP_RTI'
                error('rti_phase is only supported for SQP_RTI');
            end
            if self.outputs.KKT_residuals && strcmp(solver_options.nlp_solver_type, 'SQP_RTI')
                warning(['KKT_residuals now computes the residuals of the output iterate in SQP_RTI, ', ...
                    'this leads to increased computation time, turn off this port if it is not needed. ', ...
                    'See https://github.com/acados/acados/pull/1346.']);
            end
            % validate that all inputs/outputs are 0 or 1
            self.inputs.make_consistent();
            self.outputs.make_consistent();
            if strcmp(problem_class, 'MOCP')
                input_names = AcadosOcpSimulinkOptions.nonsupported_mocp_inputs();
                for i=1:length(input_names)
                    if self.inputs.(input_names{i})
                        warning(['Simulink inputs ', input_names{i}, ' are not supported for MOCP, turning it off.']);
                        self.inputs.(input_names{i}) = 0;
                    end
                end
            end
        end
    end

    methods (Static)
        function names = nonsupported_mocp_inputs()
            names = {'y_ref', 'lg', 'ug', 'cost_W_0', 'cost_W', 'cost_W_e'};
        end
        function obj = from_struct(data)
            obj = AcadosOcpSimulinkOptions();
            fields = fieldnames(data);
            for i = 1:length(fields)
                f = fields{i};
                try
                    obj.(f) = data.(f);
                catch
                    warning(['Could not assign field ' f ' in AcadosOcpSimulinkOptions.from_struct']);
                end
            end

            if isempty(obj.customizable_inputs)
                obj.customizable_inputs = struct();
            end
        end
    end
end