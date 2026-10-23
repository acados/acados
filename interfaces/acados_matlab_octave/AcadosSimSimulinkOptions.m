%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

classdef AcadosSimSimulinkOptions < handle

    properties
        % AcadosSimSimulinkInputs and AcadosSimSimulinkOutputs port flags.
        inputs
        outputs
    end

    methods
        function obj = AcadosSimSimulinkOptions()
            obj.inputs = AcadosSimSimulinkInputs();
            obj.outputs = AcadosSimSimulinkOutputs();
        end

        function s = to_struct(self)
            s = struct();
            s.inputs = self.inputs.to_struct();
            s.outputs = self.outputs.to_struct();
            s = orderfields(s);
        end

        function make_consistent(self, solver_options)
            self.inputs.make_consistent();
            self.outputs.make_consistent();

            if (self.inputs.xdot || self.inputs.z) && ~strcmp(solver_options.integrator_type, 'IRK')
                error('Simulink inputs xdot and z are only supported for IRK.');
            end
            if self.inputs.T && strcmp(solver_options.integrator_type, 'GNSF')
                error('Simulink input T is not supported for GNSF.');
            end
            if self.outputs.z && ~solver_options.output_z
                warning('Simulink output z requires solver_options.output_z = true, turning it off.');
                self.outputs.z = 0;
            end
            if self.outputs.S_forw && ~solver_options.sens_forw
                warning('Simulink output S_forw requires solver_options.sens_forw = true, turning it off.');
                self.outputs.S_forw = 0;
            end
            if self.outputs.S_adj && ~solver_options.sens_adj
                warning('Simulink output S_adj requires solver_options.sens_adj = true, turning it off.');
                self.outputs.S_adj = 0;
            end
        end
    end

    methods (Static)
        function obj = from_struct(data)
            obj = AcadosSimSimulinkOptions();
            names = fieldnames(data);
            for i = 1:numel(names)
                try
                    if strcmp(names{i}, 'inputs')
                        obj.inputs = AcadosSimSimulinkInputs.from_struct(data.inputs);
                    elseif strcmp(names{i}, 'outputs')
                        obj.outputs = AcadosSimSimulinkOutputs.from_struct(data.outputs);
                    else
                        warning(['Could not assign field ', names{i}, ' in AcadosSimSimulinkOptions.from_struct']);
                    end
                catch
                    warning(['Could not assign field ', names{i}, ' in AcadosSimSimulinkOptions.from_struct']);
                end
            end
        end
    end
end