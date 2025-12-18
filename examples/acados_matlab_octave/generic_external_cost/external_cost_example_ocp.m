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

import casadi.*

clear all; clc;
check_acados_requirements()

model_path = fullfile(pwd,'..','pendulum_on_cart_model');
addpath(model_path)

%% model
model = get_pendulum_on_cart_model();
nx = length(model.x);
nu = length(model.u);
ny_0 = nu; % number of outputs in initial cost term
ny = nx+nu; % number of outputs in lagrange term
ny_e = nx; % number of outputs in terminal cost term
%% discretization
N = 20;
T = 1; % time horizon length
%% set up OCP
ocp = AcadosOcp();
ocp.model = model;
ocp.model.name = 'pendulum';
sym_x = model.x;
sym_u = model.u;

params = SX.sym('Qdiag', 4, 1);
ocp.model.p = params;
ocp.parameter_values = zeros(size(params));% initialize to zero, change later

ocp.solver_options.tf = T;
ocp.solver_options.N_horizon = N;
ocp.solver_options.nlp_solver_type = 'SQP'; % 'SQP_RTI'
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_cond_N = 5; % for partial condensing
ocp.solver_options.integrator_type = 'ERK';

% cost
% generic cost
W_x = diag([1e3, 1e3, 1e-2, 1e-2]);
W_u = 1e-2;
cost_expr_ext_cost_e = 0.5 * sym_x'* W_x * sym_x;
cost_expr_ext_cost = cost_expr_ext_cost_e + 0.5 * sym_u' * W_u * sym_u;
cost_expr_ext_cost_0 = 0.5 * sym_u' * W_u * sym_u;

% nonlinear least squares
cost_expr_y_0 = sym_u;
cost_W_0 = W_u;
cost_expr_y = vertcat(sym_x, sym_u);
cost_W = blkdiag(W_x, W_u);
cost_expr_y_e = sym_x;
cost_W_e = W_x;

% linear least squares
cost_Vx_0 = zeros(ny_0,nx);
cost_Vu_0 = eye(nu);
cost_y_ref_0 = zeros(ny_0, 1);

cost_Vx = [eye(nx); zeros(nu,nx)]; % state-to-output matrix in lagrange term
cost_Vu = [zeros(nx, nu); eye(nu)]; % input-to-output matrix in lagrange term
cost_y_ref = zeros(ny, 1); % output reference in lagrange term

cost_Vx_e = eye(ny_e, nx);
cost_y_ref_e = zeros(ny_e, 1);

ocp.cost.cost_type_0 = 'EXTERNAL';
ocp.cost.cost_type = 'EXTERNAL';
ocp.cost.cost_type_e = 'EXTERNAL';

generic_or_casadi = 0; % 0=generic, 1=casadi, 2=mixed
if (generic_or_casadi == 0)
    % Generic initial cost
    ocp.cost.cost_ext_fun_type_0 = 'generic';
    ocp.cost.cost_source_ext_cost_0 = 'generic_ext_cost.c';
    ocp.cost.cost_function_ext_cost_0 =  'ext_cost';
    % Generic stage cost
    ocp.cost.cost_ext_fun_type = 'generic';
    ocp.cost.cost_source_ext_cost = 'generic_ext_cost.c';
    ocp.cost.cost_function_ext_cost =  'ext_cost';
    % Generic terminal cost
    ocp.cost.cost_ext_fun_type_e = 'generic';
    ocp.cost.cost_source_ext_cost_e = 'generic_ext_cost.c';
    ocp.cost.cost_function_ext_cost_e =  'ext_costN';
elseif (generic_or_casadi == 1)
    % Casadi initial cost
    ocp.model.cost_expr_ext_cost_0 = cost_expr_ext_cost_0;
    % Casadi stage cost
    ocp.model.cost_expr_ext_cost = cost_expr_ext_cost;
    % Casadi terminal cost
    ocp.model.cost_expr_ext_cost_e = cost_expr_ext_cost_e;
elseif (generic_or_casadi == 2)
    % Casadi initial cost
    ocp.model.cost_expr_ext_cost_0 = cost_expr_ext_cost_0;
    % Generic stage cost
    ocp.cost.cost_ext_fun_type = 'generic';
    ocp.cost.cost_source_ext_cost = 'generic_ext_cost.c';
    ocp.cost.cost_function_ext_cost =  'ext_cost';
    % Casadi terminal cost
    ocp.model.cost_expr_ext_cost_e = cost_expr_ext_cost_e;
end

% constraints
ocp.constraints.constr_type = 'AUTO';
constr_expr_h_0 = sym_u;
constr_expr_h = sym_u;
x0 = [0; pi; 0; 0];
% constraints
ocp.constraints.x0 = x0;
ocp.model.con_h_expr_0 = constr_expr_h_0;
ocp.model.con_h_expr = constr_expr_h;
U_max = 80;
ocp.constraints.lh_0 = -U_max;
ocp.constraints.uh_0 = U_max;
ocp.constraints.lh = -U_max;
ocp.constraints.uh = U_max;

%% create ocp solver
ocp_solver = AcadosOcpSolver(ocp);

x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);
% diagonal matrix Q as runtime param
p = [1e3; 1e3; 1e-2; 1e-2];

%% call ocp solver
% update initial state
ocp_solver.set('constr_x0', x0);

% set trajectory initialization
ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);
ocp_solver.set('init_pi', zeros(nx, N));

% change values for specific shooting node using:
%   ocp_solver.set('field', value, optional: stage_index)
ocp_solver.set('constr_lbx', x0, 0);

ocp_solver.set('p', p);

% solve
ocp_solver.solve();
% get solution
utraj = ocp_solver.get('u');
xtraj = ocp_solver.get('x');

status = ocp_solver.get('status'); % 0 - success
ocp_solver.print('stat')

%% Plots
figure; hold on;
States = {'p', 'theta', 'v', 'dtheta'};
for i=1:length(States)
    subplot( length(States), 1, i);
    plot(xtraj(i,:)); grid on;
    legend(States{i});
end

figure
stairs([utraj'; utraj(end)])
grid on
