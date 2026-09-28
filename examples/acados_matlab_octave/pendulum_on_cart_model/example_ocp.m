%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

import casadi.*

check_acados_requirements()


%% solver settings
N = 20; % number of discretization steps
T = 1; % [s] prediction horizon length
x0 = [0; pi; 0; 0]; % initial state

%% model dynamics
model = get_pendulum_on_cart_model();
nx = length(model.x); % state size
nu = length(model.u); % input size

%% OCP formulation object
ocp = AcadosOcp();
ocp.model = model;

%% cost in nonlinear least squares form
W_x = diag([1e3, 1e3, 1e-2, 1e-2]);
W_u = 1e-2;

% initial cost term
ny_0 = nu;
ocp.cost.cost_type_0 = 'NONLINEAR_LS';
ocp.cost.W_0 = W_u;
ocp.cost.yref_0 = zeros(ny_0, 1);
ocp.model.cost_y_expr_0 = model.u;

% path cost term
ny = nx + nu;
ocp.cost.cost_type = 'NONLINEAR_LS';
ocp.cost.W = blkdiag(W_x, W_u);
ocp.cost.yref = zeros(ny, 1);
ocp.model.cost_y_expr = vertcat(model.x, model.u);

% terminal cost term
ny_e = nx;
ocp.cost.cost_type_e = 'NONLINEAR_LS';
ocp.model.cost_y_expr_e = model.x;
ocp.cost.yref_e = zeros(ny_e, 1);
ocp.cost.W_e = W_x;

%% define constraints
% only bound on u on initial stage and path
ocp.model.con_h_expr = model.u;
ocp.model.con_h_expr_0 = model.u;

U_max = 80;
ocp.constraints.lh = -U_max;
ocp.constraints.lh_0 = -U_max;
ocp.constraints.uh = U_max;
ocp.constraints.uh_0 = U_max;
ocp.constraints.x0 = x0;

% define solver options
ocp.solver_options.N_horizon = N;
ocp.solver_options.tf = T;
ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.integrator_type = 'ERK';
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_mu0 = 1e3;
ocp.solver_options.qp_solver_cond_N = 5;
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
ocp.solver_options.globalization = 'MERIT_BACKTRACKING';
ocp.code_gen_options.ext_fun_compile_flags = '-O2';
ocp.code_gen_options.code_export_directory = 'c_generated_code'; % needed only for code reuse test
ocp.code_gen_options.json_file = 'pendulum_ocp.json'; % needed only for code reuse test

% create solver
ocp_solver = AcadosOcpSolver(ocp);

% solver initial guess
x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);

%% call ocp solver
% set trajectory initialization
ocp_solver.set('init_x', x_traj_init); % states
ocp_solver.set('init_u', u_traj_init); % inputs
ocp_solver.set('init_pi', zeros(nx, N)); % multipliers for dynamics equality constraints

% solve
ocp_solver.solve();
% get solution
utraj = ocp_solver.get('u');
xtraj = ocp_solver.get('x');

status = ocp_solver.get('status'); % 0 - success
ocp_solver.print('stat')

%% plots
ts = linspace(0, T, N+1);
figure; hold on;
states = {'p', 'theta', 'v', 'dtheta'};
for i=1:length(states)
    subplot(length(states), 1, i);
    plot(ts, xtraj(i,:)); grid on;
    ylabel(states{i});
    xlabel('t [s]')
end

figure
stairs(ts, [utraj'; utraj(end)])
ylabel('F [N]')
xlabel('t [s]')
grid on
