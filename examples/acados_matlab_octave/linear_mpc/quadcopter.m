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

% this function implements the quadcopter model from https://osqp.org/docs/examples/mpc.html
% more infomation can be found here: https://github.com/orgs/osqp/discussions/558

function model = quadcopter(h)
% input: sampling time
import casadi.*

% system matrices
Ad = [1       0       0   0   0   0   0.1     0       0    0       0       0;
      0       1       0   0   0   0   0       0.1     0    0       0       0;
      0       0       1   0   0   0   0       0       0.1  0       0       0;
      0.0488  0       0   1   0   0   0.0016  0       0    0.0992  0       0;
      0      -0.0488  0   0   1   0   0      -0.0016  0    0       0.0992  0;
      0       0       0   0   0   1   0       0       0    0       0       0.0992;
      0       0       0   0   0   0   1       0       0    0       0       0;
      0       0       0   0   0   0   0       1       0    0       0       0;
      0       0       0   0   0   0   0       0       1    0       0       0;
      0.9734  0       0   0   0   0   0.0488  0       0    0.9846  0       0;
      0      -0.9734  0   0   0   0   0      -0.0488  0    0       0.9846  0;
      0       0       0   0   0   0   0       0       0    0       0       0.9846];
Bd = [0      -0.0726  0       0.0726;
     -0.0726  0       0.0726  0;
     -0.0152  0.0152 -0.0152  0.0152;
      0      -0.0006 -0.0000  0.0006;
      0.0006  0      -0.0006  0;
      0.0106  0.0106  0.0106  0.0106;
      0      -1.4512  0       1.4512;
     -1.4512  0       1.4512  0;
     -0.3049  0.3049 -0.3049  0.3049;
      0      -0.0236  0       0.0236;
      0.0236  0      -0.0236  0;
      0.2107  0.2107  0.2107  0.2107];
[nx, nu] = size(Bd);  % state and input dimensions

% (unnamed) symbolic variables
sym_x = SX.sym('x',nx,1);  % state vector
sym_u = SX.sym('u',nu,1);  % input vector
xr = SX.sym('xr',nx,1);  % state reference
sym_p = xr;  % include the reference as a parameter

% discrete system dynamics
dyn_expr_phi = Ad * sym_x + Bd * sym_u;  % x+ = A*x + B*u

% cost matrices
Q = diag([0 0 10 10 10 10 0 0 0 5 5 5]);  % state cost
R = 0.1*eye(nu);  % input cost

% generic cost formulation
cost_expr_ext_cost_e = (sym_x-xr)'*Q*(sym_x-xr);  % terminal cost (only states)
cost_expr_ext_cost = cost_expr_ext_cost_e + sym_u'*R*sym_u;  % stage cost (states and inputs)
cost_expr_ext_cost = 1/h * cost_expr_ext_cost;  % scale the stage cost to match the discrete formulation
cost_expr_ext_cost_0 = 1/h * sym_u'*R*sym_u;  % penalize only the inputs in the first stage
% more info on discrete cost scaling: 
% https://docs.acados.org/python_interface/index.html#acados_template.acados_ocp_cost.AcadosOcpCost

% linear least-squares cost formulation (alternative)
% initial cost
ny_0 = nu;
Vu_0 = eye(ny_0);
Vx_0 = zeros(ny_0,nx);
W_0 = 1/h * 2 * R;  % scale to match the original cost
y_ref_0 = zeros(ny_0,1);

% intermediate cost
ny = nu+nx;
Vu = [eye(nu); zeros(nx,nu)];
Vx = [zeros(nu, nx); eye(nx)];
W = 1/h * 2 * blkdiag(R,Q);  % scale to match the original cost
y_ref = zeros(ny, 1);

% terminal cost
ny_e = nx;
Vx_e = eye(ny_e, nx);
W_e = 2 * Q;  % terminal cost is not scaled with h later
y_ref_e = zeros(ny_e, 1);

% input constraints
u0 = 10.5916;  % steady-state input
Jbu = eye(nu);  % all inputs are constrained
lbu = [9.6; 9.6; 9.6; 9.6] - u0;  % input lower bounds
ubu = [13; 13; 13; 13] - u0;  % input upper bounds

% state constraints on the first, second and sixth state
% (not applied to the initial state)
Jbx = zeros(3,nx);
Jbx(1,1) = 1;
Jbx(2,2) = 1;
Jbx(3,6) = 1;
infty = 1e6;  % to approximate one-sided constraints
lbx = [-pi/6; -pi/6; -1];  % state lower bounds
ubx = [ pi/6;  pi/6; infty];  % state upper bounds

% use the same constraints on the terminal state
Jbx_e = Jbx;
lbx_e = lbx;
ubx_e = ubx;

% populate the output structure
model.Ad = Ad;
model.Bd = Bd;
model.nx = nx;
model.nu = nu;

model.sym_x = sym_x;
model.sym_u = sym_u;
model.sym_p = sym_p;

model.dyn_expr_phi = dyn_expr_phi;

model.cost_expr_ext_cost_0 = cost_expr_ext_cost_0;
model.cost_expr_ext_cost = cost_expr_ext_cost;
model.cost_expr_ext_cost_e = cost_expr_ext_cost_e;

model.Vu_0 = Vu_0;
model.Vx_0 = Vx_0;
model.W_0 = W_0;
model.y_ref_0 = y_ref_0;

model.Vu = Vu;
model.Vx = Vx;
model.W = W;
model.y_ref = y_ref;

model.Vx_e = Vx_e;
model.W_e = W_e;
model.y_ref_e = y_ref_e;

model.Jbu = Jbu;
model.lbu = lbu;
model.ubu = ubu;
model.Jbx = Jbx;
model.lbx = lbx;
model.ubx = ubx;
model.Jbx_e = Jbx_e;
model.lbx_e = lbx_e;
model.ubx_e = ubx_e;
end
