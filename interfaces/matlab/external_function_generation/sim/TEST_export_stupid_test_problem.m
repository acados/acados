function [ model ] = TEST_export_stupid_test_problem()
%% this function generates an implicit ODE / index-1 DAE model,
% which consists of a CasADi expression f_impl_expr
% that depends on the symbolic CasADi variables x, xdot, u, z,
% and a model name, which will be used as a prefix for generated C
% functions later on;

%% set up f_impl (index 1 DAE)
import casadi.*

model_name_prefix = 'stupid_test';

% NOTE: this model has no physical meaning whatsoever;
% It was just made up to test stuff..

x  = SX.sym('x',5,1);     % Differential States
xdot = SX.sym('xdot', size(x));
u = SX.sym('u',3,1);
z = SX.sym('z');
p = SX.sym('parameters',2,1);

f_impl_expr = vertcat(x(1)^2 + sin(xdot(3)),...
        x(2) * 42 - pi * u(1) + xdot(2) - 1337, ...
        xdot(5) - u(3)^2 + x(4), ...
        x(2) + u(3) - xdot(1),...
        xdot(4) - u(1)  - x(2) + p(2),...
        x(2) + sqrt(x(1)) - exp(z(1)) );

model.f_impl_expr = f_impl_expr;
model.x = x;
model.xdot = xdot;
model.u = u;
model.z = z;
model.p = p;
model.name = model_name_prefix;

end