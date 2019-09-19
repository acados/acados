classdef acados_constraint < handle
    properties
        expr
        x
        u
        z
        name
        nc
    end
    methods
        function obj = acados_constraint()
        obj.expr = [];
        obj.x = [];
        obj.u = [];
        obj.z = [];
        obj.nc = [];
        obj.name = [];
        end
    end
end
    