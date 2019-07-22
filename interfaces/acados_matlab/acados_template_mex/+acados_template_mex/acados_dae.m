classdef acados_dae < handle
    properties
        f_impl_expr
        f_expl_expr
        x
        xdot
        u
        z
        name
        p
    end
    methods
        function obj = acados_dae()
        obj.f_impl_expr = [];
        obj.f_expl_expr = [];
        obj.x = [];
        obj.xdot = [];
        obj.u = [];
        obj.z = [];
        obj.name = [];
        obj.p = [];
        end
    end
end
    