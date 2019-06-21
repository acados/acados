function [ model ] = TEST_export_linear_test_problem()
%% this function generates an implicit ODE / index-1 DAE model,
% which consists of a CasADi expression f_impl_expr
% that depends on the symbolic CasADi variables x, xdot, u, z,
% and a model name, which will be used as a prefix for generated C
% functions later on;

%% set up f_impl, continuous implicit model
import casadi.*

model_name_prefix = 'linear_test';

% NOTE: this model has no physical meaning whatsoever;
% It was just made up to test stuff..
nx = 19;
nu = 2;
nz = 0;
np = 0;

x  = SX.sym('x',nx,1);     % Differential States
xdot = SX.sym('xdot', size(x));
u = SX.sym('u',nu,1);
z = SX.sym('z', nz);
p = SX.sym('parameters',np,1);

A = rand(nx,nx);
B = rand(nx,nu);

f_impl_expr = xdot - A * x - B * u;

model.f_impl_expr = f_impl_expr;
model.x = x;
model.xdot = xdot;
model.u = u;
model.z = z;
model.p = p;
model.name = model_name_prefix;

end