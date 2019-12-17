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
        import casadi.*
        obj.f_impl_expr = [];
        obj.f_expl_expr = [];
        obj.x = [];
        obj.xdot = SX.sym('xdot', 0, 0);
        obj.u = SX.sym('u', 0, 0);
        obj.z = SX.sym('z', 0, 0);
        obj.name = [];
        obj.p = SX.sym('p', 0, 0);
        end
    end
end
