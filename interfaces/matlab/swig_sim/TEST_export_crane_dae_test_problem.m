function [ model ] = TEST_export_crane_dae_test_problem()
%% this function generates an implicit ODE / index-1 DAE model,
% which consists of a CasADi expression f_impl_expr
% that depends on the symbolic CasADi variables x, xdot, u, z,
% and a model name, which will be used as a prefix for generated C
% functions later on;

import casadi.*

% this should end with '_' to get nice C function names
model_name_prefix = 'crane_dae';

%% Parameters (taken from Riens ACADO model)
tau1 = 0.012790605943772;   a1   = 0.047418203070092;
tau2 = 0.024695192379264;   a2   = 0.034087337273386;
g = 9.81;

%% Set up States & Controls
xC = SX.sym('xC');     %States
vC = SX.sym('vC');
xL = SX.sym('xL');     
vL = SX.sym('vL');
uC = SX.sym('uC');
uL = SX.sym('uL');
theta = SX.sym('theta');
omega = SX.sym('omega'); 
uCR = SX.sym('uCR');  % Controls
uLR = SX.sym('uLR');
q = SX.sym('q'); % a quadrature state

x = vertcat(xC, vC, xL, vL, uC, uL, theta, omega, q);
u = vertcat(uCR, uLR);

z = SX.sym('z',2); % define an algebraic state;

p = SX.sym('p',0,0);% define parameters p

%% explicit ODE formulation
f_expl = vertcat(vC, ...
                  - 1/tau1 * (vC - a1 * uC), ...
                  vL,...
                  - 1/tau2 * (vL - a2 * uL), ...
                  uCR,...
                  uLR,...
                  omega, ...
                  - (a1 * uCR * cos(theta) + g* sin(theta) + 2*vL*omega) / xL, ...
                  uCR^2 + xL^2); % dynamics of quadrature state x2;




%% set up equivalent implicit model
xdot = SX.sym('xdot',size(x));

f_impl_expr = [(f_expl - xdot );
            z - [((theta^2)/8 + xL);
      cos(omega + 0.1) + (xdot(3) - uCR*vL)^2 - z(2)]];

model.f_impl_expr = f_impl_expr;
model.x = x;
model.p = p;
model.xdot = xdot;
model.u = u;
model.z = z;
model.name = model_name_prefix;

end
