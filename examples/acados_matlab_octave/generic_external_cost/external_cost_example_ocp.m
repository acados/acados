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



% NOTE: `acados` currently supports both an old MATLAB/Octave interface (< v0.4.0)
% as well as a new interface (>= v0.4.0).

% THIS EXAMPLE still uses the OLD interface. If you are new to `acados` please start
% with the examples that have been ported to the new interface already.
% see https://github.com/acados/acados/issues/1196#issuecomment-2311822122)

import casadi.*

clear all; clc;
check_acados_requirements()

model_path = fullfile(pwd,'..','pendulum_on_cart_model');
addpath(model_path)

%% discretization
N = 20;
T = 1; % time horizon length
x0 = [0; pi; 0; 0];

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

%% runtime parameters
params = SX.sym('Qdiag', 4, 1);

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
ocp_model.set('sym_p', params);

% cost
ocp_model.set('cost_type_0', 'ext_cost');
ocp_model.set('cost_type', 'ext_cost');
ocp_model.set('cost_type_e', 'ext_cost');

generic_or_casadi = 0; % 0=generic, 1=casadi, 2=mixed
if (generic_or_casadi == 0)
    % Generic initial cost
    ocp_model.set('cost_ext_fun_type_0', 'generic');
    ocp_model.set('cost_source_ext_cost_0', 'generic_ext_cost.c');
    ocp_model.set('cost_function_ext_cost_0', 'ext_cost');
    % Generic stage cost
    ocp_model.set('cost_ext_fun_type', 'generic');
    ocp_model.set('cost_source_ext_cost', 'generic_ext_cost.c');
    ocp_model.set('cost_function_ext_cost', 'ext_cost');
    % Generic terminal cost
    ocp_model.set('cost_ext_fun_type_e', 'generic');
    ocp_model.set('cost_source_ext_cost_e', 'generic_ext_cost.c');
    ocp_model.set('cost_function_ext_cost_e', 'ext_costN');
elseif (generic_or_casadi == 1)
    % Casadi initial cost
    ocp_model.set('cost_ext_fun_type_0', 'casadi');
    ocp_model.set('cost_expr_ext_cost_0', model.cost_expr_ext_cost_0);
    % Casadi stage cost
    ocp_model.set('cost_ext_fun_type', 'casadi');
    ocp_model.set('cost_expr_ext_cost', model.cost_expr_ext_cost);
    % Casadi terminal cost
    ocp_model.set('cost_ext_fun_type_e', 'casadi');
    ocp_model.set('cost_expr_ext_cost_e', model.cost_expr_ext_cost_e);
elseif (generic_or_casadi == 2)
    % Casadi initial cost
    ocp_model.set('cost_ext_fun_type_0', 'casadi');
    ocp_model.set('cost_expr_ext_cost_0', model.cost_expr_ext_cost_0);
    % Generic stage cost
    ocp_model.set('cost_ext_fun_type', 'generic');
    ocp_model.set('cost_source_ext_cost', 'generic_ext_cost.c');
    ocp_model.set('cost_function_ext_cost', 'ext_cost');
    % Casadi terminal cost
    ocp_model.set('cost_ext_fun_type_e', 'casadi');
    ocp_model.set('cost_expr_ext_cost_e', model.cost_expr_ext_cost_e);
end

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
ocp_model.set('constr_expr_h_0', model.constr_expr_h_0);
ocp_model.set('constr_expr_h', model.constr_expr_h);
U_max = 80;
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
ocp_opts.set('parameter_values', zeros(size(params))); % initialize to zero, change later

%% create ocp solver
ocp_solver = acados_ocp(ocp_model, ocp_opts);

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
