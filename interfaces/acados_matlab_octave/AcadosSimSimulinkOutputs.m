%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


classdef AcadosSimSimulinkOutputs < handle

    properties
        % State at interval end, algebraic state, sensitivities, and total time.
        x
        z
        S_forw
        S_adj
        time_tot
    end

    methods
        function obj = AcadosSimSimulinkOutputs()
            obj.x = 1;
            obj.z = 0;
            obj.S_forw = 0;
            obj.S_adj = 0;
            obj.time_tot = 0;
        end

        function s = to_struct(self)
            if exist('properties')
                publicProperties = eval('properties(self)');
            else
                publicProperties = fieldnames(self);
            end
            s = struct();
            for fi = 1:numel(publicProperties)
                property_name = publicProperties{fi};
                s.(property_name) = self.(property_name);
            end
        end

        function make_consistent(self)
            s = self.to_struct();
            names = fieldnames(s);
            for i = 1:numel(names)
                value = s.(names{i});
                if ~(isequal(value, 0) || isequal(value, 1))
                    error(['AcadosSimSimulinkOutputs.', names{i}, ' must be 0 or 1, got ', num2str(value)]);
                end
            end
        end
    end

    methods (Static)
        function obj = from_struct(data)
            obj = AcadosSimSimulinkOutputs();
            names = fieldnames(data);
            for i = 1:numel(names)
                try
                    obj.(names{i}) = data.(names{i});
                catch
                    warning(['Could not assign field ', names{i}, ' in AcadosSimSimulinkOutputs.from_struct']);
                end
            end
        end
    end
end