%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


% TODO this class is not used yet!
classdef AcadosSimDims < handle
    properties

        nx     % number of states
        nu     % number of inputs
        nz     % number of algebraic variables
        np     % number of parameters
        np_global % number of global parameters, always zero for sim
    end

    methods
        function obj = AcadosSimDims()

            obj.nx = [];
            obj.nu = [];
            obj.nz = 0;
            obj.np = 0;
            obj.np_global = 0;
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
end

