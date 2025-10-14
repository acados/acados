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

%% Test of CONVEX_OVER_NONLINEAR cost type in MATLAB interface

import casadi.*

check_acados_requirements()

model_name = 'conl_cost_test';

%% Model
nx = 2;
nu = 1;

x = SX.sym('x', nx);
u = SX.sym('u', nu);

% Simple double integrator dynamics
f_expl = vertcat(x(2), u);

model = AcadosModel();
model.name = model_name;
model.x = x;
model.u = u;
model.f_expl_expr = f_expl;

%% OCP
ocp = AcadosOcp();
ocp.model = model;

% Horizon
N = 20;
T = 1.0;
ocp.solver_options.tf = T;
ocp.solver_options.N_horizon = N;

% Cost (CONVEX_OVER_NONLINEAR)
ny = nx + nu;
ny_e = nx;

% Cost type
ocp.cost.cost_type = 'CONVEX_OVER_NONLINEAR';
ocp.cost.cost_type_e = 'CONVEX_OVER_NONLINEAR';

% Cost expressions
ocp.model.cost_y_expr = vertcat(x, u);
ocp.model.cost_y_expr_e = x;

% Define residual variables
r = SX.sym('r', ny);
r_e = SX.sym('r_e', ny_e);
ocp.model.cost_r_in_psi_expr = r;
ocp.model.cost_r_in_psi_expr_e = r_e;

% Define outer convex function (quadratic for this test)
W = diag([1.0, 1.0, 0.01]);
W_e = diag([1.0, 1.0]);
ocp.model.cost_psi_expr = 0.5 * r.' * W * r;
ocp.model.cost_psi_expr_e = 0.5 * r_e.' * W_e * r_e;

% Reference
ocp.cost.yref = zeros(ny, 1);
ocp.cost.yref_e = zeros(ny_e, 1);

% Initial state
x0 = [1.0; 0.0];
ocp.constraints.x0 = x0;

% Input bounds
ocp.constraints.lbu = -2.0;
ocp.constraints.ubu = 2.0;
ocp.constraints.idxbu = 0;

% Solver options
ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_cond_N = N;
ocp.solver_options.nlp_solver_max_iter = 100;
ocp.solver_options.nlp_solver_tol_stat = 1e-6;

%% Create solver
solver = AcadosOcpSolver(ocp);

%% Solve
solver.solve();
status = solver.get('status');

if status ~= 0
    error(['Solver returned status ', num2str(status)]);
end

%% Check solution
x_sol = zeros(nx, N+1);
u_sol = zeros(nu, N);

for i = 0:N
    x_sol(:, i+1) = solver.get('x', i);
    if i < N
        u_sol(:, i+1) = solver.get('u', i);
    end
end

fprintf('\nCONVEX_OVER_NONLINEAR cost test completed successfully!\n');
fprintf('Initial state: [%.3f, %.3f]\n', x_sol(1, 1), x_sol(2, 1));
fprintf('Final state: [%.3f, %.3f]\n', x_sol(1, end), x_sol(2, end));
fprintf('Max control input: %.3f\n', max(abs(u_sol(:))));

% Check that final state is close to origin
if norm(x_sol(:, end)) > 1e-2
    error('Final state is not close to origin!');
end

fprintf('Test passed!\n');
