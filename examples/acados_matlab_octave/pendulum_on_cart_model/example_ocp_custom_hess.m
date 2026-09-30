%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

clear all; clc;
import casadi.*

% model_path = fullfile(pwd,'..','pendulum_on_cart_model');
% addpath(model_path)
% check_acados_requirements()

%% discretization
N = 40;
T = 2; % time horizon length
x0 = [0; pi; 0; 0];
xf = [0; 0; 0; 0];

nlp_solver = 'sqp'; % sqp, sqp_rti
qp_solver = 'partial_condensing_hpipm';
    % full_condensing_hpipm, partial_condensing_hpipm, full_condensing_qpoases
qp_solver_cond_N = 5; % for partial condensing
% integrator type
sim_method = 'erk'; % erk, irk, irk_gnsf

%% model dynamics
model = get_pendulum_on_cart_model();
nx = length(model.x);
nu = length(model.u);
model_name = 'pendulum';

%% OCP formulation
ocp = AcadosOcp();
ocp.model = model;
ocp.model.name = model_name;

% cost
W_u = 1e-3;
theta = model.x(2);
model.cost_expr_ext_cost = tanh(theta)^2 + .5 * (model.x(1)^2 + W_u* model.u^2);
model.cost_expr_ext_cost_e = tanh(theta)^2 + .5 * model.x(1)^2;

custom_hess_u = W_u;
% J is jacobian of inner (linear function);

J = horzcat(SX.eye(2), SX(2,2));
% diagonal matrix with second order terms of outer loss function.
D = SX.sym('D', Sparsity.diag(2));
D(1, 1) = 1;
[hess_tan, grad_tan] = hessian( tanh(theta)^2, theta);
D(2, 2) = if_else(theta == 0, hess_tan, grad_tan / theta);

custom_hess_x = J' * D * J;
if is_octave()
    error("This example does not work in Octave, somehow blkdiag doesn't work with symbolic matrices. If you happen to know how to fix this, please open a pull request.");
end
cost_expr_ext_cost_custom_hess = blkdiag(custom_hess_u, custom_hess_x);
cost_expr_ext_cost_custom_hess_e = custom_hess_x;


model.cost_expr_ext_cost_custom_hess = cost_expr_ext_cost_custom_hess;
model.cost_expr_ext_cost_custom_hess_e = cost_expr_ext_cost_custom_hess_e;
ocp.cost.cost_type = 'EXTERNAL';
ocp.cost.cost_type_e = 'EXTERNAL';

% dynamics
if (strcmp(sim_method, 'erk'))
    ocp.model.f_expl_expr = model.f_expl_expr;
    ocp.solver_options.integrator_type = 'ERK';
else % irk irk_gnsf
    ocp.model.f_impl_expr = model.f_impl_expr;
    ocp.solver_options.integrator_type = 'IRK';
end

% constraints
ocp.constraints.constr_type_0 = 'AUTO';
ocp.constraints.constr_type = 'AUTO';
ocp.model.con_h_expr_0 = model.u;
ocp.model.con_h_expr = model.u;
U_max = 35;
ocp.constraints.lh_0 = -U_max;
ocp.constraints.uh_0 = U_max;
ocp.constraints.lh = -U_max;
ocp.constraints.uh = U_max;
ocp.constraints.x0 = x0;

ocp.solver_options.N_horizon = N;
ocp.solver_options.tf = T;
ocp.solver_options.nlp_solver_type = upper(nlp_solver);
ocp.solver_options.hessian_approx = 'EXACT';
ocp.solver_options.integrator_type = upper(sim_method);
ocp.solver_options.qp_solver = upper(qp_solver);
ocp.solver_options.qp_solver_cond_N = qp_solver_cond_N;
ocp.solver_options.globalization = 'MERIT_BACKTRACKING';
ocp.solver_options.nlp_solver_max_iter = 500;

ocp_solver = AcadosOcpSolver(ocp);

x_traj_init = zeros(nx, N+1);

taus = linspace(0,1, N+1);
for i=1:N+1
    x_traj_init(:, 1) = x0*(1-taus(i)) + xf*taus(i);
end
u_traj_init = zeros(nu, N);

%% call ocp solver
% update initial state
ocp_solver.set('constr_x0', x0);

% set trajectory initialization
ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);
% ocp_solver.set('init_pi', zeros(nx, N))

% change values for specific shooting node using:
%   ocp_solver.set('field', value, optional: stage_index)
% ocp_solver.set('constr_lbx', x0, 0)

% solve
ocp_solver.solve();
% get solution
utraj = ocp_solver.get('u');
xtraj = ocp_solver.get('x');

status = ocp_solver.get('status'); % 0 - success
ocp_solver.print('stat')

%% Plots
ts = linspace(0, T, N+1);
figure; hold on;
States = {'p', 'theta', 'v', 'dtheta'};
for i=1:length(States)
    subplot(length(States), 1, i);
    plot(ts, xtraj(i,:)); grid on;
    ylabel(States{i});
    xlabel('t [s]')
end

figure
stairs(ts, [utraj'; utraj(end)])
ylabel('F [N]')
xlabel('t [s]')
grid on
