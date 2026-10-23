%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

classdef AcadosSimSimulinkInputs < handle

    properties
        % Initial state, control, and parameter ports. x/u/p default on.
        x
        u
        p
        % IRK-only initial guesses for xdot and algebraic state z.
        xdot
        z
        % Adjoint seed (nx), simulation interval T (except GNSF), and start time t0.
        seed_adj
        T
        t0
    end

    methods
        function obj = AcadosSimSimulinkInputs()
            obj.x = 1;
            obj.u = 1;
            obj.p = 1;
            obj.xdot = 0;
            obj.z = 0;
            obj.seed_adj = 0;
            obj.T = 0;
            obj.t0 = 0;
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
                    error(['AcadosSimSimulinkInputs.', names{i}, ' must be 0 or 1, got ', num2str(value)]);
                end
            end
        end
    end

    methods (Static)
        function obj = from_struct(data)
            obj = AcadosSimSimulinkInputs();
            names = fieldnames(data);
            for i = 1:numel(names)
                try
                    obj.(names{i}) = data.(names{i});
                catch
                    warning(['Could not assign field ', names{i}, ' in AcadosSimSimulinkInputs.from_struct']);
                end
            end
        end
    end
end