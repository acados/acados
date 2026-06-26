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

clear all
addpath('.');

% Check env.sh
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
    warning('env.sh does not seem to be sourced. If the build fails, run: source env.sh');
end

%% Load model
raw = ocp_model_wind_turbine_nx6();

nx = raw.nx;
nu = raw.nu;
np = raw.np;

model = AcadosModel();
model.name = 'wind_turbine_nx6';
model.x = raw.sym_x;
model.u = raw.sym_u;
model.xdot = raw.sym_xdot;
model.p = raw.sym_p;
model.f_expl_expr = raw.expr_f_expl;
model.f_impl_expr = raw.expr_f_impl;

%% Set up OCP object
ocp = AcadosOcp();
ocp.model = model;

% Horizon & shooting nodes
N = 40;
Ts = 0.2; % Sampling time
T = N * Ts; % Horizon length [s]

%% Cost (Linear least squares)
ny = 4; % Outputs in Lagrange term
ny_e = 2; % Outputs in Mayer term

ocp.cost.cost_type = 'LINEAR_LS';
ocp.cost.cost_type_e = 'LINEAR_LS';
ocp.cost.cost_type_0 = 'LINEAR_LS';

% State-to-output matrix (Lagrange)
Vx = zeros(ny, nx);
Vx(1, 1) = 1.0;
Vx(2, 5) = 1.0;

% Input-to-output matrix (Lagrange)
Vu = zeros(ny, nu);
Vu(3, 1) = 1.0;
Vu(4, 2) = 1.0;

% State-to-output matrix (Mayer)
Vx_e = zeros(ny_e, nx);
Vx_e(1, 1) = 1.0;
Vx_e(2, 5) = 1.0;

% Weight matrix (Lagrange)
W = zeros(ny, ny);
W(1, 1) = 1.5114;
W(2, 1) = -0.0649;
W(1, 2) = -0.0649;
W(2, 2) = 0.0180;
W(3, 3) = 0.01;
W(4, 4) = 0.001;

% Weight matrix (Mayer)
W_e = zeros(ny_e, ny_e);
W_e(1, 1) = 1.5114;
W_e(2, 1) = -0.0649;
W_e(1, 2) = -0.0649;
W_e(2, 2) = 0.0180;

% Initial-node cost
ocp.cost.Vx_0 = Vx;
ocp.cost.Vu_0 = Vu;
ocp.cost.W_0 = W;
ocp.cost.yref_0 = zeros(ny, 1);

% Intermediate cost
ocp.cost.Vx = Vx;
ocp.cost.Vu = Vu;
ocp.cost.W = W;
ocp.cost.yref = zeros(ny, 1);

% Terminal cost
ocp.cost.Vx_e = Vx_e;
ocp.cost.W_e = W_e;
ocp.cost.yref_e = zeros(ny_e, 1);

% Slack penalties
Z = 1e2;
z = 0e2;
Z_e = 1e2;
z_e = 0e2;

% Intermediate soft nonlinear constraint
ocp.cost.Zl = Z;
ocp.cost.Zu = Z;
ocp.cost.zl = z;
ocp.cost.zu = z;

% Terminal soft nonlinear constraint
ocp.cost.Zl_e = Z_e;
ocp.cost.Zu_e = Z_e;
ocp.cost.zl_e = z_e;
ocp.cost.zu_e = z_e;

use_soft_h0_constraint = 1;
if use_soft_h0_constraint
    % Initial-node soft nonlinear constraint
    ocp.cost.Zl_0 = Z;
    ocp.cost.Zu_0 = Z;
    ocp.cost.zl_0 = z;
    ocp.cost.zu_0 = z;
end

%% Constraints
% Bound constants
dbeta_min = -8.0;
dbeta_max = 8.0;
dM_gen_min = -1.0;
dM_gen_max = 1.0;
OmegaR_min = 6.0 / 60 * 2 * pi;
OmegaR_max = 13.0 / 60 * 2 * pi;
beta_min = 0.0;
beta_max = 35.0;
M_gen_min = 0.0;
M_gen_max = 5.0;
Pel_min = 0.0;
Pel_max = 5.0;

% State bounds (0-based indices)
ocp.constraints.idxbx = [0; 6; 7];
ocp.constraints.lbx = [OmegaR_min; beta_min; M_gen_min];
ocp.constraints.ubx = [OmegaR_max; beta_max; M_gen_max];

% Input bounds
ocp.constraints.idxbu = [0; 1];
ocp.constraints.lbu = [dbeta_min; dM_gen_min];
ocp.constraints.ubu = [dbeta_max; dM_gen_max];

% Nonlinear power constraint
ocp.model.con_h_expr = raw.expr_h;
ocp.constraints.lh = Pel_min;
ocp.constraints.uh = Pel_max;

ocp.model.con_h_expr_e = raw.expr_h_e;
ocp.constraints.lh_e = Pel_min;
ocp.constraints.uh_e = Pel_max;

if use_soft_h0_constraint
    ocp.model.con_h_expr_0 = raw.expr_h;
    ocp.constraints.lh_0 = Pel_min;
    ocp.constraints.uh_0 = Pel_max;
end

% Soft nonlinear constraints
ocp.constraints.idxsh = 0;
ocp.constraints.idxsh_e = 0;
if use_soft_h0_constraint
    ocp.constraints.idxsh_0 = 0;
end

% Initial state
ocp.constraints.x0 = zeros(nx, 1);

% Parameter default value
ocp.parameter_values = zeros(np, 1);

%% Create solver
ocp.solver_options.N_horizon = N;
ocp.solver_options.tf = T;

% Solver options
ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
ocp.solver_options.regularize_method = 'NO_REGULARIZE';
ocp.solver_options.nlp_solver_ext_qp_res = 1;
ocp.solver_options.nlp_solver_max_iter = 200;

ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_cond_N = 5;
ocp.solver_options.qp_solver_cond_ric_alg = 0;
ocp.solver_options.qp_solver_ric_alg = 0;
ocp.solver_options.qp_solver_warm_start = 0;
ocp.solver_options.qp_solver_iter_max = 50;

% Integrator
ocp.solver_options.integrator_type = 'IRK';
ocp.solver_options.sim_method_num_stages = 4;
ocp.solver_options.sim_method_num_steps = 1;

ocp_solver = AcadosOcpSolver(ocp);

%% References / Initialization
compute_setup;

% Trajectory initialization
x_traj_init = repmat(x0_ref, 1, N + 1);
u_traj_init = repmat(u0_ref, 1, N);
ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);

% Set x0
ocp_solver.set('constr_x0', x0_ref);

% Set time-varying parameter per stage
for jj = 0:N
    ocp_solver.set('p', wind0_ref(:, jj + 1), jj);
end

% Set references per stage
for jj = 0:N - 1
    ocp_solver.set('cost_y_ref', y_ref(:, jj + 1), jj);
end
ocp_solver.set('cost_y_ref', y_ref(1:ny_e, N + 1), N);

%% Solve
disp('before solve')
tic;
ocp_solver.solve();
time_ext = toc;
disp('after solve')

% Get solution
u = ocp_solver.get('u');
x = ocp_solver.get('x');

electrical_power = 0.944 * 97 / 100 * x(1, :) .* x(6, :);

%% Diagnostics
status = ocp_solver.get('status');
sqp_iter = ocp_solver.get('sqp_iter');
time_tot = ocp_solver.get('time_tot');
time_lin = ocp_solver.get('time_lin');
time_reg = ocp_solver.get('time_reg');
time_qp_sol = ocp_solver.get('time_qp_sol');

fprintf(['\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] ', ...
    '(time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n'], ...
    status, sqp_iter, time_ext * 1e3, time_tot * 1e3, ...
    time_lin * 1e3, time_qp_sol * 1e3, time_reg * 1e3);

ocp_solver.print('stat');

%% Figures
figure;

subplot(3, 1, 1);
plot(0:N, x);
xlim([0, N]);
ylabel('states');

subplot(3, 1, 2);
plot(0:N - 1, u);
xlim([0, N]);
ylabel('inputs');

subplot(3, 1, 3);
plot(0:N, electrical_power);
hold on;
plot([0, N], [Pel_max, Pel_max]);
hold off;
xlim([0, N]);
ylim([4.5, 6.0]);
ylabel('electrical power');

stat = ocp_solver.get('stat');
figure;
plot(0:size(stat, 1) - 1, log10(stat(:, 2)), 'r-x');
hold on;
plot(0:size(stat, 1) - 1, log10(stat(:, 3)), 'b-x');
plot(0:size(stat, 1) - 1, log10(stat(:, 4)), 'g-x');
plot(0:size(stat, 1) - 1, log10(stat(:, 5)), 'k-x');
hold off;
xlabel('iter');
ylabel('res');

if status == 0
    fprintf('\nsuccess!\n\n');
else
    fprintf('\nsolution failed!\n\n');
end

if is_octave()
    waitforbuttonpress;
end