function [ model ] = TEST_simple_model()
%% this function generates an explicit ODE test model,
% represented as a Matlab struct "model".
% It consists of a CasADi expression f_expl_expr
% that depends on the symbolic CasADi variables x, xdot, u (also part of
% "model"), a model name, which will be used as a prefix for the name of 
% generated C functions to use the model with acados.

import casadi.*

model_name = 'crane_dae';

%% Parameters
tau1 = 0.012790605943772;   a1   = 0.047418203070092;
tau2 = 0.024695192379264;   a2   = 0.034087337273386;
g = 9.81;

%% Set up States & Controls
xC = MX.sym('xC');     %States
vC = MX.sym('vC');
xL = MX.sym('xL');     
vL = MX.sym('vL');
uC = MX.sym('uC');
uL = MX.sym('uL');
theta = MX.sym('theta');
omega = MX.sym('omega');
q = MX.sym('q');

uCR = MX.sym('uCR');  % Controls
uLR = MX.sym('uLR');

x = vertcat(xC, vC, xL, vL, uC, uL, theta, omega, q);
u = vertcat(uCR, uLR);

xdot = MX.sym('xdot',size(x)); %state derivatives

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
model.f_expl_expr = f_expl;
model.f_impl_expr = f_expl - xdot;
model.x = x;
model.xdot = xdot;
model.u = u;
model.name = model_name;

end
