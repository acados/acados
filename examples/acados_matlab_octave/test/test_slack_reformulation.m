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



function test_slack_reformulation()

import casadi.*

%% settings
N       = 20;
Tf      = 1.0;
Nsim    = 20;
Fmax    = 80;
vmax    = 5;
tol_closed_loop = 5e-5;
plot_result = false;

% qp_solvers = {'PARTIAL_CONDENSING_HPIPM', ...
%               'FULL_CONDENSING_QPOASES', ...
%               'FULL_CONDENSING_HPIPM', ...
%               'FULL_CONDENSING_DAQP', ...
%               'PARTIAL_CONDENSING_OSQP'};

qp_solvers = {'PARTIAL_CONDENSING_HPIPM', ...
              'FULL_CONDENSING_QPOASES', ...
              'FULL_CONDENSING_HPIPM', ...
              'FULL_CONDENSING_DAQP'};

settings = struct('N', N, 'Tf', Tf, 'Nsim', Nsim, 'Fmax', Fmax, 'vmax', vmax, ...
                  'plot_result', plot_result);

%% 1) slack formulations (bx vs h) for all QP solvers
types = {'bx', 'h'};
results = cell(numel(types), numel(qp_solvers));
labels = {};
for i = 1:numel(types)
    for j = 1:numel(qp_solvers)
        lbl = sprintf('%s | %s', types{i}, qp_solvers{j});
        res = run_closed_loop(types{i}, qp_solvers{j}, false, settings);
        results{i, j} = res;
        labels{end+1} = lbl; %#ok<SAGROW>
    end
end

ref = results{1, 1};
for k = 1:numel(labels)
    res = results{ceil(k / numel(qp_solvers)), mod(k - 1, numel(qp_solvers)) + 1};
    compare_results(ref, res, labels{k}, tol_closed_loop);
end
fprintf('soft constraint example: SUCCESS, equivalent formulations agree up to %.2e\n', tol_closed_loop);

%% 2) explicit quadratic penalty vs. slack formulation (quadratic-only penalty)
ref_h   = run_closed_loop('h',         'PARTIAL_CONDENSING_HPIPM', true, settings);
res_pen = run_closed_loop('h_penalty', 'PARTIAL_CONDENSING_HPIPM', true, settings);
compare_results(ref_h, res_pen, 'penalty formulation', tol_closed_loop);
fprintf('soft constraint example: SUCCESS, penalty formulation matches slack formulation\n');

end


%% ======================================================================
%% local functions
%% ======================================================================

function model = export_pendulum_ode_model()
    model = AcadosModel();
    model.name = 'pendulum';

    % parameters
    M = 1.0;    % mass of cart [kg]
    m = 0.1;    % mass of ball [kg]
    g = 9.81;   % gravity [m/s^2]
    l = 0.8;    % length of rod [m]

    % states / controls
    x1     = casadi.SX.sym('x1');
    theta  = casadi.SX.sym('theta');
    v1     = casadi.SX.sym('v1');
    dtheta = casadi.SX.sym('dtheta');
    x      = vertcat(x1, theta, v1, dtheta);
    u      = casadi.SX.sym('F');
    xdot   = casadi.SX.sym('xdot', 4, 1);

    denom = M + m - m*cos(theta)^2;
    f_expl = vertcat( ...
        v1, ...
        dtheta, ...
        (-m*l*sin(theta)*dtheta^2 + m*g*cos(theta)*sin(theta) + u) / denom, ...
        (-m*l*cos(theta)*sin(theta)*dtheta^2 + u*cos(theta) + (M+m)*g*sin(theta)) / (l*denom));

    model.x = x;
    model.xdot = xdot;
    model.u = u;
    model.f_expl_expr = f_expl;
    model.f_impl_expr = xdot - f_expl;
end


function res = run_closed_loop(soft_constr_type, qp_solver, quadratic_penalty_only, s)
    ocp = AcadosOcp();
    model = export_pendulum_ode_model();
    ocp.model = model;
    ocp.name = [model.name soft_constr_type mat2str(quadratic_penalty_only)];

    nx = length(model.x);
    nu = length(model.u);
    ny = nx + nu;

    ocp.solver_options.N_horizon = s.N;
    ocp.solver_options.tf = s.Tf;

    % ---- cost
    Q_mat = 2*diag([1e3, 1e3, 1e-2, 1e-2]);
    R_mat = 2*diag(1e-2);
    W     = blkdiag(Q_mat, R_mat);
    Zl = 10; Zu = 10;
    % Zl = 0;
    % Zu = 0;

    x0 = [0.0; pi; 0.0; 0.0];
    ocp.constraints.x0 = x0;

    if strcmp(soft_constr_type, 'h_penalty')
        % explicit penalty in an external cost, no slack variables
        v   = model.x(3);
        y   = vertcat(model.x, model.u);
        y_e = model.x;
        penalty = 0.5*Zl*fmax(0, -s.vmax - v)^2 + 0.5*Zu*fmax(0, v - s.vmax)^2;

        ocp.cost.cost_type   = 'EXTERNAL';
        ocp.cost.cost_type_e = 'EXTERNAL';
        ocp.model.cost_expr_ext_cost   = 0.5*(y'*W*y) + penalty;
        ocp.model.cost_expr_ext_cost_e = 0.5*(y_e'*Q_mat*y_e);
    else
        ocp.cost.cost_type   = 'LINEAR_LS';
        ocp.cost.cost_type_e = 'LINEAR_LS';
        ocp.cost.W   = W;
        ocp.cost.W_e = Q_mat;
        ocp.cost.Vx  = [eye(nx); zeros(nu, nx)];
        ocp.cost.Vu  = [zeros(nx, nu); eye(nu)];
        ocp.cost.Vx_e = eye(nx);
        ocp.cost.yref   = zeros(ny, 1);
        ocp.cost.yref_e = zeros(nx, 1);

        lin_weight = 50;
        if quadratic_penalty_only
            lin_weight = 0;
        end
        ocp.cost.zl = lin_weight;
        ocp.cost.zu = lin_weight;
        ocp.cost.Zl = Zl;
        ocp.cost.Zu = Zu;
    end

    % ---- constraints
    ocp.constraints.idxbu = 0;          % 0-based indices in the template interface
    ocp.constraints.lbu = -s.Fmax;
    ocp.constraints.ubu = +s.Fmax;

    switch soft_constr_type
        case 'bx'
            ocp.constraints.idxbx  = 2;   % v is x(3) -> index 2 (0-based)
            ocp.constraints.lbx    = -s.vmax;
            ocp.constraints.ubx    = +s.vmax;
            ocp.constraints.idxsbx = 0;
        case 'h'
            ocp.model.con_h_expr = model.x(3);
            ocp.constraints.lh = -s.vmax;
            ocp.constraints.uh = +s.vmax;
            ocp.constraints.idxsh = 0;
        case 'h_penalty'
            % nothing else to do
        otherwise
            error("soft_constr_type must be 'bx', 'h' or 'h_penalty', got %s", soft_constr_type);
    end

    % ---- solver options
    ocp.solver_options.qp_solver = qp_solver;
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
    ocp.solver_options.integrator_type = 'ERK';
    ocp.solver_options.nlp_solver_type = 'SQP';
    ocp.solver_options.nlp_solver_tol_stat = 1e-7;
    ocp.solver_options.nlp_solver_tol_eq = 1e-7;
    ocp.solver_options.nlp_solver_tol_ineq = 1e-7;
    ocp.solver_options.nlp_solver_tol_comp = 1e-7;
    ocp.solver_options.nlp_solver_ext_qp_res = 1;
    ocp.solver_options.qp_solver_warm_start = 0;
    ocp.solver_options.qp_solver_iter_max = 10000;

    ocp_solver = AcadosOcpSolver(ocp);

    % ---- closed loop
    simX = zeros(s.Nsim+1, nx);
    simU = zeros(s.Nsim, nu);
    sqp_iter = zeros(s.Nsim, 1);
    x = x0;

    for i = 1:s.Nsim
        simX(i, :) = x';

        ocp_solver.set('constr_x0', x);
        ocp_solver.solve();
        status = ocp_solver.get('status');
        if status ~= 0
            ocp_solver.print('stat');
            error('OCP solver returned status %d at closed-loop step %d.', status, i);
        end

        sqp_iter(i) = ocp_solver.get('sqp_iter');
        simU(i, :)  = ocp_solver.get('u', 0)';

        % sim_solver.set('x', x);
        % sim_solver.set('u', simU(i, :)');
        % sim_solver.solve();
        % x = sim_solver.get('xn');
        x = ocp_solver.get('x', 1);
    end
    simX(s.Nsim+1, :) = x';

    % slack values at stage 1 (only for slack-based formulations)
    if any(strcmp(soft_constr_type, {'bx', 'h'}))
        fprintf('sl %g, su %g\n', ocp_solver.get('sl', 1), ocp_solver.get('su', 1));
    end

    if s.plot_result
        t = linspace(0, s.Nsim*s.Tf/s.N, s.Nsim+1);
        figure; subplot(3,1,1); stairs(t(1:end-1), simU); ylabel('F [N]'); grid on;
        subplot(3,1,2); plot(t, simX(:,1)); ylabel('x [m]'); grid on;
        subplot(3,1,3); plot(t, simX(:,2)); ylabel('theta [rad]'); xlabel('t [s]'); grid on;
    end

    fprintf('\nsoft constraint example: formulation %s with %s ran successfully.\n', ...
            soft_constr_type, qp_solver);
    fprintf('took SQP iterations: %s\n', mat2str(sqp_iter'));

    res = struct('simX', simX, 'simU', simU, 'sqp_iter', sqp_iter);
    clear ocp_solver
end


function compare_results(ref, res, label, tol)
    error_x  = norm(ref.simX - res.simX);
    error_u  = norm(ref.simU - res.simU);
    error_xu = max(error_x, error_u);
    fprintf('soft constraint example: %s deviates from reference by %.3e\n', label, error_xu);

    if error_xu > tol
        error('%s: solutions should match up to %.2e, got error_x %g, error_u %g.', ...
              label, tol, error_x, error_u);
    end
    if any(ref.sqp_iter ~= res.sqp_iter)
        error('%s: all formulations should take the same number of SQP iterations.', label);
    end
end