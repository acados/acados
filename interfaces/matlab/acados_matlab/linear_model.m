%% this function generates an explicit ODE test model,
% represented as a Matlab struct "model".
% It consists of a CasADi expression f_expl_expr
% that depends on the symbolic CasADi variables x, xdot, u (also part of
% "model"), a model name, which will be used as a prefix for the name of 
% generated C functions to use the model with acados.

import casadi.*

nx = 4;
nu = 0;

x = MX.sym('x', nx, 1); % states
u = MX.sym('u', nu, nu); % controls
xdot = MX.sym('xdot',size(x)); %state derivatives

expr_expl = [1.0; -1.0; 0.5; 0.1].*x;
expr_impl = expr_expl - xdot;

