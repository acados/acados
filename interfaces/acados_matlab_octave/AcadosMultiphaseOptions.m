%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.



classdef AcadosMultiphaseOptions < handle
    properties
        integrator_type
        collocation_type
        cost_discretization
    end
    methods
        function obj = AcadosMultiphaseOptions()
            obj.integrator_type = {};
            obj.collocation_type = {};
            obj.cost_discretization = {};
        end

        function make_consistent(self, opts, n_phases)
            % check if all fields are lists of length n_phases
            INTEGRATOR_TYPE_VALUES = {'ERK', 'ERK_WITH_COST', 'IRK', 'GNSF', 'DISCRETE', 'LIFTED_IRK'};
            COLLOCATION_TYPE_VALUES = {'GAUSS_RADAU_IIA', 'GAUSS_LEGENDRE', 'EXPLICIT_RUNGE_KUTTA'};
            COST_DISCRETIZATION_VALUES = {'EULER', 'INTEGRATOR'};

            prop_names = {'integrator_type', 'collocation_type', 'cost_discretization'};
            for prop = prop_names
                prop = prop{1};
                if ~iscell(self.(prop))
                    error('AcadosMultiphaseOptions.%s must be a cell array, got %s.', prop, class(self.(prop)));
                end
                if isempty(self.(prop))
                    % non varying field, use value from ocp opts
                    for i = 1:n_phases
                        self.(prop){i} = opts.(prop);
                    end
                elseif length(self.(prop)) ~= n_phases
                    error('AcadosMultiphaseOptions.%s must be a cell array of length n_phases, got %d.', prop, length(self.(prop)));
                end
                for i = 1:n_phases
                    if ~any(strcmp(self.(prop){i}, eval([upper(prop), '_VALUES'])))
                        error('AcadosMultiphaseOptions.%s{%d} must be one of %s, got %s.', prop, i, strjoin(eval([upper(prop), '_VALUES']), ', '), self.(prop){i});
                    end
                end
            end
        end
        function s = to_struct(self)
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
    methods (Static)
        function obj = from_struct(s)
            % Create AcadosMultiphaseOptions from a struct (e.g. decoded from JSON).
            obj = AcadosMultiphaseOptions();
            fields = fieldnames(s);
            for i = 1:length(fields)
                f = fields{i};
                % direct assignment for simple fields
                try
                    obj.(f) = s.(f);
                catch
                    % ignore unknown fields
                    warning(['Could not assign field ' f ' in AcadosMultiphaseOptions.from_struct']);
                end
            end
        end
    end
end
