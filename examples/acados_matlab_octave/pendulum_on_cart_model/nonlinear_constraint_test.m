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


clear all; clc;

check_acados_requirements()

%% OCP DESCRIPTION
ocp = AcadosOcp();

%% SOLVER OPTIONS

%% discretization
N = 10; % 40
T = 2.0; % time horizon length
h = T/N;

% uniform time grid
time_steps = T/N * ones(N,1);

shooting_nodes = zeros(N+1, 1);
for i = 1:N
    shooting_nodes(i+1) = sum(time_steps(1:i));
end

ocp.solver_options.tf = T;
ocp.solver_options.N_horizon = N;
ocp.solver_options.time_steps = time_steps;
ocp.solver_options.nlp_solver_type = 'SQP'; % 'SQP_RTI'
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'; % 'EXACT'
ocp.solver_options.regularize_method = 'CONVEXIFY';
% NO_REGULARIZE, PROJECT, PROOJECT_REDUC_HESS, MIRROR, CONVEXIFY
ocp.solver_options.nlp_solver_max_iter = 50;
ocp.solver_options.nlp_solver_tol_stat = 1e-8;
ocp.solver_options.nlp_solver_tol_eq = 1e-8;
ocp.solver_options.nlp_solver_tol_ineq = 1e-8;
ocp.solver_options.nlp_solver_tol_comp = 1e-8;
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_cond_N = 5; % for partial condensing
ocp.solver_options.qp_solver_cond_ric_alg = 0;
ocp.solver_options.qp_solver_ric_alg = 0;
ocp.solver_options.qp_solver_warm_start = 1; % 0: cold, 1: warm, 2: hot
ocp.solver_options.qp_solver_iter_max = 1000; % default is 50; OSQP needs a lot sometimes.
ocp.solver_options.qp_solver_mu0 = 1e4;
ocp.solver_options.print_level = 1;
ocp.solver_options.store_iterates = true;
ocp.solver_options.log_dual_step_norm = true;
ocp.solver_options.log_primal_step_norm = true;

% can vary for integrators
sim_method_num_stages = 1 * ones(N,1);
sim_method_num_stages(3:end) = 2;
ocp.solver_options.sim_method_num_stages = sim_method_num_stages;
ocp.solver_options.sim_method_num_steps = ones(N,1);

% integrator type
integrator = 1;
switch integrator
case 1
    ocp.solver_options.integrator_type = 'ERK';
case 2
    ocp.solver_options.integrator_type = 'IRK';
case 3
    if ~all(time_steps == T/N)
        error('nonuniform time discretization with discrete dynamics should not be used');
    end
    ocp.solver_options.integrator_type = 'DISCRETE';
otherwise
    ocp.solver_options.integrator_type = 'GNSF';
end

%% MODEL
model = get_pendulum_on_cart_model();
ocp.model = model;

% dimensions
nx = model.x.rows();
nu = model.u.rows();


%% COST
cost_formulation = 1;
switch cost_formulation
case 1
    cost_type = 'LINEAR_LS';
case 2
    cost_type = 'EXTERNAL';
otherwise
    cost_type = 'AUTO';
end

ocp.cost.cost_type_0 = cost_type;
ocp.cost.cost_type = cost_type;
ocp.cost.cost_type_e = cost_type;

W_x = diag([1e3, 1e3, 1e-2, 1e-2]);
W_u = 1e-2;

cost_expr_ext_cost_e = 0.5 * model.x'* W_x * model.x;
cost_expr_ext_cost = cost_expr_ext_cost_e + 0.5 * model.u' * W_u * model.u;
cost_expr_ext_cost_0 = 0.5 * model.u' * W_u * model.u;

ny_0 = nu; % number of outputs in initial cost term
Vx_0 = zeros(ny_0,nx);
Vu_0 = eye(nu);
y_ref_0 = zeros(ny_0, 1);

ny = nx+nu; % number of outputs in lagrange term
Vx = [eye(nx); zeros(nu,nx)]; % state-to-output matrix in lagrange term
Vu = [zeros(nx, nu); eye(nu)]; % input-to-output matrix in lagrange term
y_ref = zeros(ny, 1); % output reference in lagrange term

ny_e = nx; % number of outputs in terminal cost term
Vx_e = eye(ny_e, nx);
y_ref_e = zeros(ny_e, 1);

if strcmp(cost_type, 'LINEAR_LS')
    ocp.cost.Vu_0 = Vu_0;
    ocp.cost.Vx_0 = Vx_0;
    ocp.cost.W_0 = W_u;
    ocp.cost.yref_0 = y_ref_0;

    ocp.cost.Vu = Vu;
    ocp.cost.Vx = Vx;
    ocp.cost.W = blkdiag(W_x, W_u);
    ocp.cost.yref = y_ref;

    ocp.cost.Vx_e = Vx_e;
    ocp.cost.W_e = W_x;
    ocp.cost.yref_e = y_ref_e;
else % EXTERNAL, AUTO
    ocp.model.cost_expr_ext_cost_0 = cost_expr_ext_cost_0;
    ocp.model.cost_expr_ext_cost = cost_expr_ext_cost;
    ocp.model.cost_expr_ext_cost_e = cost_expr_ext_cost_e;
end

%% CONSTRAINTS
constraint_formulation_nonlinear = 0;
lbu = -80*ones(nu, 1);
ubu =  80*ones(nu, 1);


% bound on u
ocp.constraints.idxbu = [0];
ocp.constraints.lbu = lbu;
ocp.constraints.ubu = ubu;

ocp.constraints.idxbx = [0];
ocp.constraints.lbx = -3;
ocp.constraints.ubx = 3;

% nonlinear constraint
infty = get_acados_infty();
ocp.model.con_h_expr = model.u*(model.x(1)+3.5);
ocp.constraints.lh = -infty;
% ocp.constraints.lh = lbu;
ocp.constraints.uh = ubu;


% initial state
x0 = [0; pi; 0; 0];
ocp.constraints.x0 = x0;

%% SOLVER
ocp_solver = AcadosOcpSolver(ocp);

%% INITIALIZATION
x_traj_init = zeros(nx, N+1);
x_traj_init(2, :) = linspace(pi, 0, N+1); % initialize theta

u_traj_init = zeros(nu, N);

%% SOLVE
% prepare evaluation
n_executions = 1;
time_tot = zeros(n_executions,1);
time_lin = zeros(n_executions,1);
time_reg = zeros(n_executions,1);
time_qp_sol = zeros(n_executions,1);

%% call ocp solver in loop
for i=1:n_executions
    % initial state
    ocp_solver.set('constr_x0', x0);

    % set trajectory initialization
    ocp_solver.set('init_x', x_traj_init);
    ocp_solver.set('init_u', u_traj_init);
    ocp_solver.set('init_pi', zeros(nx, N)); % multipliers for dynamics equality constraints

    % solve
    ocp_solver.solve();
    % get solution
    utraj = ocp_solver.get('u');
    xtraj = ocp_solver.get('x');

    % evaluation
    status = ocp_solver.get('status');
    sqp_iter = ocp_solver.get('sqp_iter');
    time_tot(i) = ocp_solver.get('time_tot');
    time_lin(i) = ocp_solver.get('time_lin');
    time_reg(i) = ocp_solver.get('time_reg');
    time_qp_sol(i) = ocp_solver.get('time_qp_sol');

    if i == 1 || i == n_executions
        ocp_solver.print('stat');
    end
end

% test constraint evaluation
stat_mat = ocp_solver.get('stat');
res_ineq = stat_mat(:, 4);
flag_test = 0;
for iteration = 1:3
    % load iterate
    iterate = ocp_solver.get_iterate(iteration);
    ocp_solver.load_iterate_from_obj(iterate);

    % expected infeasibility
    ineq_res_as_iter = res_ineq(iteration+1);

    % evaluate constraints
    ineq_fun = ocp_solver.evaluate_constraints_and_get_violation();
    violations = zeros(N+1, 1);
    for i=1:length(ineq_fun)
        violations(i) = max([ineq_fun{i}; 0]);
    end
    for i=2:N
        if ineq_fun{i}(3) ~= 0.0 % 3 is lh index
            error('inequality value corresponding to masked constraint should be 0.0.');
        end
    end
    [max_violation, index] = max(violations);
    if abs(max_violation - ineq_res_as_iter) > 1e-6
        error('inequality constraint violation does not match expected value');
    else
        fprintf('inequality constraint violation matches expected value: %f\n', max_violation);
    end
    violation_idx = ocp_solver.get_constraint_indices_with_violation(max_violation);
    if max_violation ~= 0.0
        if ~isequal(size(violation_idx), [1, 2])
            error('expected violation index to be of size [1, 2], got %d, %d', size(violation_idx));
        end
        if ineq_fun{violation_idx(1)+1}(violation_idx(2)+1) ~= max_violation
            error('inequality constraint violation does not match expected value');
        else
            fprintf('max. constraint violation index correctly identified as: %d %d\n', violation_idx(1), violation_idx(2));
        end
        flag_test = 1;
    end
end

if ~flag_test
    error('constraint violation index test was not done');
end

% test res_stat getter
res_stat_all = ocp_solver.get('res_stat_all');
res_stat_norm = stat_mat(end, 2);
max_res_stat = zeros(1, length(res_stat_all));
if length(res_stat_all) ~= N+1
    error('length of res_stat_all does not match N+1');
end
flag_test = 0;
for i=1:length(res_stat_all)
    max_res_stat(i) = max(res_stat_all{i});
    if max_res_stat(i) > res_stat_norm
        error('max res_stat is larger than res_stat_norm');
    elseif max_res_stat(i) == res_stat_norm
        flag_test = 1;
    end
end

if ~flag_test
    error('did not find max res_stat equal to res_stat_norm');
else
    disp('Found max res_stat equal to res_stat_norm');
end

% test res_eq getter
res_eq_all = ocp_solver.get('res_eq_all');
res_eq_norm = stat_mat(end, 3);
max_res_eq = zeros(1, length(res_eq_all));
if length(res_eq_all) ~= N
    error('length of res_eq_all does not match N');
end
flag_test = 0;
for i=1:length(res_eq_all)
    max_res_eq(i) = max(res_eq_all{i});
    if max_res_eq(i) > res_eq_norm
        error('max res_eq is larger than res_eq_norm');
    elseif max_res_eq(i) == res_eq_norm
        flag_test = 1;
    end
end

if ~flag_test
    error('did not find max res_eq equal to res_eq_norm');
else
    disp('Found max res_eq equal to res_eq_norm');
end



%% Plot trajectories
if 0
    figure; hold on;
    States = {'p', 'theta', 'v', 'dtheta'};
    for i=1:length(States)
        subplot(length(States), 1, i);
        plot(shooting_nodes, xtraj(i,:)); grid on;
        ylabel(States{i});
        xlabel('t [s]')
    end

    figure
    stairs(shooting_nodes, [utraj'; utraj(end)])

    ylabel('F [N]')
    xlabel('t [s]')
    grid on
    if is_octave()
        waitforbuttonpress;
    end
end
