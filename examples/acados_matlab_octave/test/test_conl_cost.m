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
%% This test compares CONVEX_OVER_NONLINEAR with NONLINEAR_LS cost types
%% to verify they produce the same solution for a quadratic cost

import casadi.*

check_acados_requirements()

%% Common problem setup
nx = 2;
nu = 1;
N = 20;
T = 1.0;
x0 = [1.0; 0.0];
ny = nx + nu;
ny_e = nx;

% Cost matrices (quadratic cost)
W = diag([1.0, 1.0, 0.01]);
W_e = diag([1.0, 1.0]);

% References
yref = zeros(ny, 1);
yref_e = zeros(ny_e, 1);

%% Test 1: CONVEX_OVER_NONLINEAR formulation
fprintf('\n=== Test 1: CONVEX_OVER_NONLINEAR cost ===\n');

x = SX.sym('x', nx);
u = SX.sym('u', nu);
f_expl = vertcat(x(2), u);

model_conl = AcadosModel();
model_conl.name = 'conl_cost_test';
model_conl.x = x;
model_conl.u = u;
model_conl.f_expl_expr = f_expl;

ocp_conl = AcadosOcp();
ocp_conl.model = model_conl;
ocp_conl.solver_options.tf = T;
ocp_conl.solver_options.N_horizon = N;

% CONVEX_OVER_NONLINEAR cost setup
ocp_conl.cost.cost_type = 'CONVEX_OVER_NONLINEAR';
ocp_conl.cost.cost_type_e = 'CONVEX_OVER_NONLINEAR';

ocp_conl.model.cost_y_expr = vertcat(x, u);
ocp_conl.model.cost_y_expr_e = x;

r = SX.sym('r', ny);
r_e = SX.sym('r_e', ny_e);
ocp_conl.model.cost_r_in_psi_expr = r;
ocp_conl.model.cost_r_in_psi_expr_e = r_e;
ocp_conl.model.cost_psi_expr = 0.5 * r.' * W * r;
ocp_conl.model.cost_psi_expr_e = 0.5 * r_e.' * W_e * r_e;

ocp_conl.cost.yref = yref;
ocp_conl.cost.yref_e = yref_e;

ocp_conl.constraints.x0 = x0;
ocp_conl.constraints.lbu = -2.0;
ocp_conl.constraints.ubu = 2.0;
ocp_conl.constraints.idxbu = 0;

ocp_conl.solver_options.nlp_solver_type = 'SQP';
ocp_conl.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp_conl.solver_options.qp_solver_cond_N = N;
ocp_conl.solver_options.nlp_solver_max_iter = 100;
ocp_conl.solver_options.nlp_solver_tol_stat = 1e-6;

solver_conl = AcadosOcpSolver(ocp_conl);
solver_conl.solve();
status_conl = solver_conl.get('status');

if status_conl ~= 0
    error(['CONVEX_OVER_NONLINEAR solver returned status ', num2str(status_conl)]);
end

% Extract solution
x_sol_conl = zeros(nx, N+1);
u_sol_conl = zeros(nu, N);
for i = 0:N
    x_sol_conl(:, i+1) = solver_conl.get('x', i);
    if i < N
        u_sol_conl(:, i+1) = solver_conl.get('u', i);
    end
end

fprintf('CONVEX_OVER_NONLINEAR solution:\n');
fprintf('  Final state: [%.6f, %.6f]\n', x_sol_conl(1, end), x_sol_conl(2, end));
fprintf('  Max control: %.6f\n', max(abs(u_sol_conl(:))));

%% Test 2: NONLINEAR_LS formulation (equivalent quadratic cost)
fprintf('\n=== Test 2: NONLINEAR_LS cost (for comparison) ===\n');

x = SX.sym('x', nx);
u = SX.sym('u', nu);
f_expl = vertcat(x(2), u);

model_nls = AcadosModel();
model_nls.name = 'nls_cost_test';
model_nls.x = x;
model_nls.u = u;
model_nls.f_expl_expr = f_expl;

ocp_nls = AcadosOcp();
ocp_nls.model = model_nls;
ocp_nls.solver_options.tf = T;
ocp_nls.solver_options.N_horizon = N;

% NONLINEAR_LS cost setup (equivalent to CONVEX_OVER_NONLINEAR with quadratic psi)
ocp_nls.cost.cost_type = 'NONLINEAR_LS';
ocp_nls.cost.cost_type_e = 'NONLINEAR_LS';

ocp_nls.model.cost_y_expr = vertcat(x, u);
ocp_nls.model.cost_y_expr_e = x;

ocp_nls.cost.W = W;
ocp_nls.cost.W_e = W_e;
ocp_nls.cost.yref = yref;
ocp_nls.cost.yref_e = yref_e;

ocp_nls.constraints.x0 = x0;
ocp_nls.constraints.lbu = -2.0;
ocp_nls.constraints.ubu = 2.0;
ocp_nls.constraints.idxbu = 0;

ocp_nls.solver_options.nlp_solver_type = 'SQP';
ocp_nls.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp_nls.solver_options.qp_solver_cond_N = N;
ocp_nls.solver_options.nlp_solver_max_iter = 100;
ocp_nls.solver_options.nlp_solver_tol_stat = 1e-6;

solver_nls = AcadosOcpSolver(ocp_nls);
solver_nls.solve();
status_nls = solver_nls.get('status');

if status_nls ~= 0
    error(['NONLINEAR_LS solver returned status ', num2str(status_nls)]);
end

% Extract solution
x_sol_nls = zeros(nx, N+1);
u_sol_nls = zeros(nu, N);
for i = 0:N
    x_sol_nls(:, i+1) = solver_nls.get('x', i);
    if i < N
        u_sol_nls(:, i+1) = solver_nls.get('u', i);
    end
end

fprintf('NONLINEAR_LS solution:\n');
fprintf('  Final state: [%.6f, %.6f]\n', x_sol_nls(1, end), x_sol_nls(2, end));
fprintf('  Max control: %.6f\n', max(abs(u_sol_nls(:))));

%% Compare solutions
fprintf('\n=== Comparing solutions ===\n');
tol = 1e-5;

x_diff = max(max(abs(x_sol_conl - x_sol_nls)));
u_diff = max(max(abs(u_sol_conl - u_sol_nls)));

fprintf('Max state difference: %.2e\n', x_diff);
fprintf('Max control difference: %.2e\n', u_diff);

if x_diff > tol
    error(['State solutions differ by more than tolerance! Max diff: ', num2str(x_diff)]);
end

if u_diff > tol
    error(['Control solutions differ by more than tolerance! Max diff: ', num2str(u_diff)]);
end

% Check that final state is close to origin
if norm(x_sol_conl(:, end)) > 1e-2
    error('Final state is not close to origin!');
end

fprintf('\n=== All tests passed! ===\n');
fprintf('CONVEX_OVER_NONLINEAR and NONLINEAR_LS produce equivalent solutions.\n');
