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
        function obj = AcadosOcpSimulinkOptions()
            obj.outputs = AcadosOcpSimulinkOutputs();
            obj.inputs = AcadosOcpSimulinkInputs();
            obj.samplingtime = 't0';
            obj.show_port_info = 1;
            obj.zoro_iterations = 1;
            obj.generate_simulink_block = 1;
            obj.customizable_inputs = struct();
        end

        function add_customizable_input(self, input_name, input_spec)
            self.customizable_inputs.(input_name) = input_spec;
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
    end

    methods (Static)
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