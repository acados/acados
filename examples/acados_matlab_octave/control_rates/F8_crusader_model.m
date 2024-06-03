%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;

%

function model = F8_crusader_model()
import casadi.*

%% System dimension
nx = 4;
nu = 1;

%% System parameters
a1 = -0.877; a2 = -0.088; a3 = 0.47; a4 = -0.019; a5 = 3.846; a6 = -0.215;
a7 = 0.28; a8 = 0.47; a9 = 0.63;
c1 = -4.208; c2 = -0.396; c3 = -0.47; c4 = -3.564; c5 = -20.967; c6 = 6.265;
c7 = 46; c8 = 61.4;

%% Symbolic variables
% states
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
x4 = SX.sym('x4');  % extra state for the original control input
sym_x = vertcat(x1, x2, x3, x4);

% controls
u = SX.sym('u');  % input rate becomes the input
sym_u = u;

% state derivatives
sym_xdot = SX.sym('xdot', nx, 1); 

%% The system dynamics
% Explicit Dynamics
expr_f_expl = vertcat(a1*x1+x3+a2*x1*x3+a3*x1.^2+a4*x2.^2-x1.^2*x3+a5*x1.^3+a6*x4+a7*x1.^2*x4+a8*x1*x4.^2+a9*x4.^3, ...
                      x3, ...
                      c1*x1+c2*x3+c3*x1.^2+c4*x1.^3+c5*x4+c6*x1.^2*x4+c7*x1*x4.^2+c8*x4.^3, ...
                      u);  % dynamics of the added state

% Implicit Dynamics
expr_f_impl = expr_f_expl - sym_xdot;

%% fill structure
model.nx = nx;
model.nu = nu;
model.sym_x = sym_x;
model.sym_u = sym_u;
model.sym_xdot = sym_xdot;
model.expr_f_expl = expr_f_expl;
model.expr_f_impl = expr_f_impl;

end