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

addpath('../pendulum_on_cart_model')

check_acados_requirements()

print_level = 1;

%% SOLVER OPTIONS
solver_options = AcadosOcpOptions();

%% discretization
N = 40;
T = 2.0; % time horizon length
h = T/N;

% nonuniform time grid
% N1 = 5;
% N2 = N - N1;
% time_steps = [( 1 * ones(N1,1)); 3 * ones(N2,1)];
% time_steps = T/sum(time_steps) * time_steps;

% uniform time grid
time_steps = T/N * ones(N,1);

shooting_nodes = zeros(N+1, 1);
for i = 1:N
    shooting_nodes(i+1) = sum(time_steps(1:i));
end

solver_options.tf = T;
solver_options.N_horizon = N;
solver_options.time_steps = time_steps;
solver_options.nlp_solver_type = 'SQP'; % 'SQP_RTI'
solver_options.hessian_approx = 'GAUSS_NEWTON'; % 'EXACT'
solver_options.regularize_method = 'CONVEXIFY';
% NO_REGULARIZE, PROJECT, PROOJECT_REDUC_HESS, MIRROR, CONVEXIFY
solver_options.nlp_solver_max_iter = 50;
solver_options.nlp_solver_tol_stat = 1e-8;
solver_options.nlp_solver_tol_eq = 1e-8;
solver_options.nlp_solver_tol_ineq = 1e-8;
solver_options.nlp_solver_tol_comp = 1e-8;
solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
% FULL_CONDENSING_HPIPM, PARTIAL_CONDENSING_HPIPM
% FULL_CONDENSING_QPOASES, PARTIAL_CONDENSING_OSQP
solver_options.qp_solver_cond_N = 5; % for partial condensing
solver_options.qp_solver_cond_ric_alg = 0;
solver_options.qp_solver_ric_alg = 0;
solver_options.qp_solver_warm_start = 1; % 0: cold, 1: warm, 2: hot
solver_options.qp_solver_iter_max = 1000; % default is 50; OSQP needs a lot sometimes.
solver_options.qp_solver_mu0 = 1e4;
solver_options.exact_hess_dyn = 1;
solver_options.exact_hess_cost = 1;
solver_options.exact_hess_constr = 1;

% can vary for integrators
sim_method_num_stages = 1 * ones(N,1);
sim_method_num_stages(3:end) = 2;
solver_options.sim_method_num_stages = sim_method_num_stages;
solver_options.sim_method_num_steps = ones(N,1);

% integrator type
integrator = 1;
switch integrator
case 1
    solver_options.integrator_type = 'ERK';
case 2
    solver_options.integrator_type = 'IRK';
case 3
    if ~all(time_steps == T/N)
        error('nonuniform time discretization with discrete dynamics should not be used');
    end
    solver_options.integrator_type = 'DISCRETE';
otherwise
    solver_options.integrator_type = 'GNSF';
end

%% MODEL
model = get_pendulum_on_cart_AcadosModel(T/N);

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

cost = AcadosOcpCost();

cost.cost_type_0 = cost_type;
cost.cost_type = cost_type;
cost.cost_type_e = cost_type;

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
    cost.Vu_0 = Vu_0;
    cost.Vx_0 = Vx_0;
    cost.W_0 = W_u;
    cost.yref_0 = y_ref_0;

    cost.Vu = Vu;
    cost.Vx = Vx;
    cost.W = blkdiag(W_x, W_u);
    cost.yref = y_ref;

    cost.Vx_e = Vx_e;
    cost.W_e = W_x;
    cost.yref_e = y_ref_e;
else % EXTERNAL, AUTO
    cost.cost_expr_ext_cost_0 = cost_expr_ext_cost_0;
    cost.cost_expr_ext_cost = cost_expr_ext_cost;
    cost.cost_expr_ext_cost_e = cost_expr_ext_cost_e;
end

%% CONSTRAINTS
constraints = AcadosOcpConstraints();
constraint_formulation_nonlinear = 0;
lbu = -80*ones(nu, 1);
ubu =  80*ones(nu, 1);

constraints.constr_type = 'AUTO';
constraints.constr_type_0 = 'AUTO';
constraints.constr_type_e = 'AUTO';

if constraint_formulation_nonlinear % formulate constraint via h
    model.con_h_expr_0 = model.u;
    constraints.lh_0 = lbu;
    constraints.uh_0 = ubu;
    constraints.con_h_expr = model.u;
    constraints.lh = lbu;
    constraints.uh = ubu;
else % formulate constraint as bound on u
    constraints.idxbu = [0];
    constraints.lbu = lbu;
    constraints.ubu = ubu;
end

% initial state
x0 = [0; pi; 0; 0];
constraints.x0 = x0;

%% OCP DESCRIPTION
ocp = AcadosOcp();
ocp.model = model;
ocp.constraints = constraints;
ocp.cost = cost;
ocp.solver_options = solver_options;

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

ocp_solver.set('print_level', print_level)

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
        ocp_solver.print('stat')
    end
end

% get slack values
for i = 0:N-1
    sl = ocp_solver.get('sl', i);
    su = ocp_solver.get('su', i);
end
sl = ocp_solver.get('sl', N);
su = ocp_solver.get('su', N);

% get cost value
cost_val_ocp = ocp_solver.get_cost();

%% get QP matrices:
% See https://docs.acados.org/problem_formulation
%        |----- dynamics -----|------ cost --------|---------------------------- constraints ------------------------|
fields = {'qp_A','qp_B','qp_b','qp_R','qp_Q','qp_r','qp_C','qp_D','qp_lg','qp_ug','qp_lbx','qp_ubx','qp_lbu','qp_ubu'};

% eiter stage-wise ...
for stage = [0,N-1]
    for k = 1:length(fields)
        field = fields{k};
        disp(strcat(field, " at stage ", num2str(stage), " = "));
        ocp_solver.get(field, stage)
    end
end

stage = N;
field = 'qp_Q';
disp(strcat(field, " at stage ", num2str(stage), " = "));
ocp_solver.get(field, stage)
field = 'qp_R';
disp(strcat(field, " at stage ", num2str(stage), " = "));
ocp_solver.get(field, stage)

... or for all stages.
qp_Q = ocp_solver.get('qp_Q');
cond_H = ocp_solver.get('qp_solver_cond_H');

disp('QP diagnostics of last QP before condensing')
result = ocp_solver.qp_diagnostics(false);
disp(['min eigenvalues of blocks are in [', num2str(min(result.min_ev)), ', ', num2str(max(result.min_ev)), ']'])
disp(['max eigenvalues of blocks are in [', num2str(min(result.max_ev)), ', ', num2str(max(result.max_ev)), ']'])
disp(['condition_number_blockwise: '])
disp(result.condition_number_blockwise)
disp(['condition_number_global: ', num2str(result.condition_number_global)])

disp('QP diagnostics of last QP after partial condensing')
result = ocp_solver.qp_diagnostics(true);
disp(['min eigenvalues of blocks are in [', num2str(min(result.min_ev)), ', ', num2str(max(result.min_ev)), ']'])
disp(['max eigenvalues of blocks are in [', num2str(min(result.max_ev)), ', ', num2str(max(result.max_ev)), ']'])
disp(['condition_number_blockwise: '])
disp(result.condition_number_blockwise)
disp(['condition_number_global: ', num2str(result.condition_number_global)])

%% Plot trajectories
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

%% plot average compuation times
% if ~is_octave()
%     time_total = sum(time_tot);
%     time_linearize = sum(time_lin);
%     time_regulariz = sum(time_reg);
%     time_qp_solution = sum(time_qp_sol);
%
%     figure;
%
%     bar_vals = 1000 * [time_linearize; time_regulariz; time_qp_solution; ...
%         time_total - time_linearize - time_regulariz - time_qp_solution] / n_executions;
%     bar([1; nan], [bar_vals, nan(size(bar_vals))]' ,'stacked')
%     legend('linearization', 'regularization', 'qp solution', 'remaining')
%     ylabel('time in [ms]')
%     title( [ strrep(cost_type, '_',' '), ' , sim: ' strrep(sim_method, '_',' '), ...
%        ';  ', strrep(qp_solver, '_', ' ')] )
% end
