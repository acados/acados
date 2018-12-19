function model = linear_model()

%% this function generates an explicit ODE test model,
% represented as a Matlab struct "model".
% It consists of a CasADi expression f_expl_expr
% that depends on the symbolic CasADi variables x, xdot, u (also part of
% "model"), a model name, which will be used as a prefix for the name of 
% generated C functions to use the model with acados.

import casadi.*

nx = 4;
nu = 0;

sym_x = MX.sym('x', nx, 1); % states
sym_xdot = MX.sym('xdot',size(sym_x)); %state derivatives

expr_f_expl = [1.0; -1.0; 0.5; 0.1].*sym_x;
expr_f_impl = expr_f_expl - sym_xdot;

model.nx = nx;
model.nu = nu;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.expr_f_expl = expr_f_expl;
model.expr_f_impl = expr_f_impl;

