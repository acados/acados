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
            INTEGRATOR_TYPE_VALUES = {'ERK', 'IRK', 'GNSF', 'DISCRETE', 'LIFTED_IRK'};
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
end
