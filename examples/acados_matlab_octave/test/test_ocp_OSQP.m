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

%% test of native matlab interface

import casadi.*
addpath('../pendulum_on_cart_model/');


%% discretization
N = 20;
T = 1; % time horizon length
x0 = [0; pi; 0; 0];

nlp_solver = 'sqp'; % sqp, sqp_rti
qp_solver = 'partial_condensing_osqp';
    % full_condensing_hpipm, partial_condensing_hpipm, full_condensing_qpoases, partial_condensing_osqp
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
ocp_model.set('constr_type', 'auto');
ocp_model.set('constr_expr_h', model.expr_h);
U_max = 80;
ocp_model.set('constr_lh', -U_max); % lower bound on h
ocp_model.set('constr_uh', U_max);  % upper bound on h

ocp_model.set('constr_x0', x0);
% ... see ocp_model.model_struct to see what other fields can be set

%% acados ocp set opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('sim_method', sim_method);
ocp_opts.set('qp_solver', qp_solver);
ocp_opts.set('qp_solver_iter_max', 2000);
ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
ocp_opts.set('ext_fun_compile_flags', '');
% ... see ocp_opts.opts_struct to see what other fields can be set

%% create ocp solver
ocp = acados_ocp(ocp_model, ocp_opts);

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
% get solutionn
utraj = ocp.get('u');
xtraj = ocp.get('x');

status = ocp.get('status'); % 0 - success
ocp.print('stat')

if status == 0
    disp('test_ocp_OSQP: success!');
else
    disp(['test_ocp_OSQP: Failed with status ', num2str(status)]);
    quit(1);
end

