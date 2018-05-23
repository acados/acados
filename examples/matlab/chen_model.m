function [ode_fun, nx, nu] = chen_model()

import casadi.*
% The following ODE model comes from Chen1998
nx = 2;
nu = 1;
x = SX.sym('x', nx);
u = SX.sym('u', nu);
mu = 0.5;
rhs = vertcat(x(2) + u*(mu + (1.-mu)*x(1)), x(1) + u*(mu - 4.*(1.-mu)*x(2)));
ode_fun = Function('chen', {x, u}, {rhs});

end
