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
x0 = [1.0; 0.5];
ny = nx + nu;
ny_e = nx;

% Cost matrices (quadratic cost)
W = diag([1.0, 1.0, 0.01]);
W_e = diag([1.0, 1.0]);

% References
yref = zeros(ny, 1);
yref_e = zeros(ny_e, 1);

%% Solve with both cost formulations
fprintf('\n=== Testing CONVEX_OVER_NONLINEAR and NONLINEAR_LS cost formulations ===\n');

cost_types = {'CONVEX_OVER_NONLINEAR', 'NONLINEAR_LS'};
x_sols = cell(2, 1);
u_sols = cell(2, 1);

for i = 1:2
    cost_type = cost_types{i};
    fprintf('\n--- Solving with %s cost ---\n', cost_type);

    % Create model
    x = SX.sym('x', nx);
    u = SX.sym('u', nu);
    f_expl = vertcat(x(2), u);

    model = AcadosModel();
    model.name = [lower(cost_type), '_test'];
    model.x = x;
    model.u = u;
    model.f_expl_expr = f_expl;

    % Create OCP
    ocp = AcadosOcp();
    ocp.model = model;
    ocp.solver_options.tf = T;
    ocp.solver_options.N_horizon = N;

    ocp.json_file = ['ocp' model.name '.json'];

    % Set cost type
    ocp.cost.cost_type = cost_type;
    ocp.cost.cost_type_e = cost_type;

    % Set cost expressions
    ocp.model.cost_y_expr = vertcat(x, u);
    ocp.model.cost_y_expr_e = x;

    ocp.cost.yref_e = yref_e;
    ocp.cost.yref = yref;

    if strcmp(cost_type, 'CONVEX_OVER_NONLINEAR')
        % CONVEX_OVER_NONLINEAR cost setup
        r = SX.sym('r', ny);
        r_e = SX.sym('r_e', ny_e);
        ocp.model.cost_r_in_psi_expr = r;
        ocp.model.cost_r_in_psi_expr_e = r_e;
        ocp.model.cost_psi_expr = 0.5 * r' * W * r;
        ocp.model.cost_psi_expr_e = 0.5 * r_e' * W_e * r_e;
    else
        % NONLINEAR_LS cost setup (equivalent to CONVEX_OVER_NONLINEAR with quadratic psi)
        ocp.cost.W = W;
        ocp.cost.W_e = W_e;
    end

    % Set constraints
    ocp.constraints.x0 = x0;
    ocp.constraints.lbu = -2.0;
    ocp.constraints.ubu = 2.0;
    ocp.constraints.idxbu = 0;

    % Set solver options
    ocp.solver_options.nlp_solver_type = 'SQP';
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
    ocp.solver_options.qp_solver_cond_N = N;
    ocp.solver_options.nlp_solver_max_iter = 100;
    ocp.solver_options.nlp_solver_tol_stat = 1e-6;

    % Create solver and solve
    solver = AcadosOcpSolver(ocp);
    solver.solve();
    status = solver.get('status');

    solver.print_statistics();

    if status ~= 0
        error(['%s solver returned status ', num2str(status)], cost_type);
    end

    % Extract solution
    x_sols{i} = solver.get('x');
    u_sols{i} = solver.get('u');
end

%% Compare solutions
fprintf('\n=== Comparing solutions ===\n');
tol = 1e-9;

x_diff = max(max(abs(x_sols{1} - x_sols{2})));
u_diff = max(max(abs(u_sols{1} - u_sols{2})));

fprintf('Max state difference: %.2e\n', x_diff);
fprintf('Max control difference: %.2e\n', u_diff);

if x_diff > tol
    error(['State solutions differ by more than tolerance! Max diff: ', num2str(x_diff)]);
end

if u_diff > tol
    error(['Control solutions differ by more than tolerance! Max diff: ', num2str(u_diff)]);
end

fprintf('\n=== All tests passed! ===\n');
fprintf('CONVEX_OVER_NONLINEAR and NONLINEAR_LS produce equivalent solutions.\n');

clear solver