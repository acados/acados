%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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

%% test of native matlab interface
clear all

model_path = fullfile(pwd,'..','pendulum_on_cart_model');
addpath(model_path)

check_acados_requirements()

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
model = pendulum_on_cart_model;
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
ocp_model.set('cost_expr_ext_cost', model.expr_ext_cost);
ocp_model.set('cost_expr_ext_cost_e', model.expr_ext_cost_e);

% dynamics
if (strcmp(sim_method, 'erk'))
    ocp_model.set('dyn_type', 'explicit');
    ocp_model.set('dyn_expr_f', model.expr_f_expl);
else % irk irk_gnsf
    ocp_model.set('dyn_type', 'implicit');
    ocp_model.set('dyn_expr_f', model.expr_f_impl);
end

% constraints
ocp_model.set('constr_type', 'auto'); % bgh
ocp_model.set('constr_expr_h', model.expr_h);
U_max = 80;
ocp_model.set('constr_lh', -U_max); % lower bound on h
ocp_model.set('constr_uh', U_max);  % upper bound on h

nbx = 1;
Jbx = zeros(nbx,nx);
Jbx(1,1) = 1;
ocp_model.set('constr_Jbx', Jbx);
ocp_model.set('constr_lbx', -1)
ocp_model.set('constr_ubx', 1)

ocp_model.set('constr_x0', x0);
% ... see ocp_model.model_struct to see what other fields can be set

%% acados ocp set opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('sim_method', sim_method);
ocp_opts.set('qp_solver', qp_solver);
ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
% ... see ocp_opts.opts_struct to see what other fields can be set

%% create ocp solver
ocp = acados_ocp(ocp_model, ocp_opts);

%% Use custom cost
generic_or_casadi = 0; % 0=generic, 1=casadi

if (generic_or_casadi == 0)
    ocp.C_ocp_ext_fun = ocp_set_ext_fun_manual_cost(ocp.C_ocp, ocp.C_ocp_ext_fun, ocp.model_struct, ocp.opts_struct, ...
        'generic', ...
        'cost/generic_ext_cost.c', 'ext_cost',...
        'cost/generic_ext_cost.c', 'ext_costN');
else
    ocp.C_ocp_ext_fun = ocp_set_ext_fun_manual_cost(ocp.C_ocp, ocp.C_ocp_ext_fun, ocp.model_struct, ocp.opts_struct, ...
        'casadi', ...
        'cost/casadi_ext_cost.c', 'pendulum_cost_ext_cost_fun_jac_hess',...
        'cost/casadi_ext_cost_e.c', 'pendulum_cost_ext_cost_e_fun_jac_hess');
end

ocp.C_ocp_ext_fun.cost_ext_cost_fun_jac_hess

% precompute
ocp_precompute(ocp.C_ocp);

%%
x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);

%% call ocp solver
% update initial state
ocp.set('constr_x0', x0);

% set trajectory initialization
ocp.set('init_x', x_traj_init);
ocp.set('init_u', u_traj_init);
ocp.set('init_pi', zeros(nx, N))

% change values for specific shooting node using:
%   ocp.set('field', value, optional: stage_index)
ocp.set('constr_lbx', x0, 0)

% solve
ocp.solve();
% get solution
utraj = ocp.get('u');
xtraj = ocp.get('x');

status = ocp.get('status'); % 0 - success
ocp.print('stat')

%% Plots
figure; hold on;
States = {'p', 'theta', 'v', 'dtheta'};
for i=1:length(States)
    subplot( length(States), 1, i);
    plot(xtraj(i,:)); grid on;
    legend(States{i});
end

figure
stairs(utraj')
grid on

%% go embedded
% to generate templated C code
% ocp.generate_c_code;
