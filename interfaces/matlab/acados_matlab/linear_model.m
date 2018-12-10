function [ model ] = linear_model(model_name)

import casadi.*

nx = 4;
nu = 0;

x = MX.sym('x', nx, 1); % states
u = MX.sym('u', nu, nu); % controls
xdot = MX.sym('xdot',size(x)); %state derivatives

f_expl = [1.0; -1.0; 0.5; 0.1].*x;

%% populate model
model.f_expl_expr = f_expl;
model.f_impl_expr = f_expl - xdot;
model.x = x;
model.xdot = xdot;
model.u = u;
model.name = model_name;
model.nx = nx;
model.nu = nu;


