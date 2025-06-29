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

% NOTE: `acados` currently supports both an old MATLAB/Octave interface (< v0.4.0)
% as well as a new interface (>= v0.4.0).

% THIS EXAMPLE still uses the OLD interface. If you are new to `acados` please start
% with the examples that have been ported to the new interface already.
% see https://github.com/acados/acados/issues/1196#issuecomment-2311822122)


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
model = pendulum_on_cart_model();
nx = model.nx;
nu = model.nu;

%% model to create the solver
ocp_model = acados_ocp_model();
model_name = 'pendulum';

%% acados ocp model
ocp_model.set('name', model_name);
ocp_model.set('T', T);

% symbolics
ocp_model.set('sym_x', model.sym_x);
ocp_model.set('sym_u', model.sym_u);
ocp_model.set('sym_xdot', model.sym_xdot);

% cost
ocp_model.set('cost_type', 'ext_cost');
ocp_model.set('cost_type_e', 'ext_cost');

W_u = 1e-3;
theta = model.sym_x(2);
model.expr_ext_cost = tanh(theta)^2 + .5 * (model.sym_x(1)^2 + W_u* model.sym_u^2);
model.expr_ext_cost_e = tanh(theta)^2 + .5 * model.sym_x(1)^2;

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


ocp_model.set('cost_expr_ext_cost', model.expr_ext_cost);
ocp_model.set('cost_expr_ext_cost_e', model.expr_ext_cost_e);
ocp_model.set('cost_expr_ext_cost_custom_hess', cost_expr_ext_cost_custom_hess);
ocp_model.set('cost_expr_ext_cost_custom_hess_e', cost_expr_ext_cost_custom_hess_e);

% dynamics
if (strcmp(sim_method, 'erk'))
    ocp_model.set('dyn_type', 'explicit');
    ocp_model.set('dyn_expr_f', model.dyn_expr_f_expl);
else % irk irk_gnsf
    ocp_model.set('dyn_type', 'implicit');
    ocp_model.set('dyn_expr_f', model.dyn_expr_f_impl);
end

% constraints
ocp_model.set('constr_type', 'auto');
ocp_model.set('constr_expr_h_0', model.constr_expr_h);
ocp_model.set('constr_expr_h', model.constr_expr_h);
U_max = 35;
ocp_model.set('constr_lh_0', -U_max); % lower bound on h
ocp_model.set('constr_uh_0', U_max);  % upper bound on h
ocp_model.set('constr_lh', -U_max);
ocp_model.set('constr_uh', U_max);

ocp_model.set('constr_x0', x0);

%% acados ocp set opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('sim_method', sim_method);
ocp_opts.set('qp_solver', qp_solver);
ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
ocp_opts.set('globalization', 'merit_backtracking');
ocp_opts.set('nlp_solver_max_iter', 500);

%% create ocp solver
ocp_solver = acados_ocp(ocp_model, ocp_opts);

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
