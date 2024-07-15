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

clear all;
check_acados_requirements()

import casadi.*

%% solver settings
N = 20; % number of discretization steps
T = 1; % [s] prediction horizon length
nlp_solver = 'sqp'; % sqp, sqp_rti
% integrator type
sim_method = 'erk'; % erk, irk, irk_gnsf

nx = 3; % state size
nu = 3; % input size
ny = nx+nu;
ny_e = nx;

x0 = ones(nx, 1); % initial state

%% acados ocp model
ocp_model = acados_ocp_model();

sym_x = SX.sym('x', nx);
sym_u = SX.sym('u', nu);


ocp_model.set('name', 'trivial');
ocp_model.set('T', T);  % prediction horizon

% symbolics
ocp_model.set('sym_x', sym_x);
ocp_model.set('sym_u', sym_u);

% cost
cost_type = 'linear_ls';
cost_type_e = 'linear_ls';
Vx = zeros(ny,nx); Vx(1:nx,:) = eye(nx);        % state-to-output matrix in lagrange term
Vu = zeros(ny,nu); Vu(nx+1:ny,:) = eye(nu);     % input-to-output matrix in lagrange term
Vx_e = zeros(ny_e,nx); Vx_e(1:nx,:) = eye(nx);  % state-to-output matrix in mayer term
W = eye(ny);
W_e = 5 * W(1:ny_e,1:ny_e);                         % cost weights in mayer term
y_ref = zeros(ny,1);                            % set intial references
y_ref_e = zeros(ny_e,1);

% cost
ocp_model.set('cost_type', cost_type);
ocp_model.set('cost_type_e', cost_type_e);
ocp_model.set('cost_Vx', Vx);
ocp_model.set('cost_Vu', Vu);
ocp_model.set('cost_Vx_e', Vx_e);
ocp_model.set('cost_W', W);
ocp_model.set('cost_W_e', W_e);
ocp_model.set('cost_y_ref', y_ref);
ocp_model.set('cost_y_ref_e', y_ref_e);

% dynamics
if (strcmp(sim_method, 'erk'))
    ocp_model.set('dyn_type', 'explicit');
    ocp_model.set('dyn_expr_f', sym_u);
end

% constraints
ocp_model.set('constr_x0', x0);  % set the initial state

%% acados ocp options
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('sim_method', sim_method);


%% Simulink opts
simulink_opts = get_acados_simulink_opts;

simulink_opts.inputs.y_ref = 0;
simulink_opts.inputs.y_ref_0 = 0;
simulink_opts.inputs.y_ref_e = 0;

simulink_opts.inputs.cost_W = 0;
simulink_opts.inputs.cost_W_0 = 0;
simulink_opts.inputs.cost_W_e = 0;

simulink_opts.inputs.x_init = 1;
simulink_opts.inputs.u_init = 1;
simulink_opts.inputs.pi_init = 1;

simulink_opts.inputs.reset_solver = 0;

simulink_opts.outputs.u0 = 0;
simulink_opts.outputs.utraj = 1;
simulink_opts.outputs.xtraj = 1;

simulink_opts.outputs.solver_status = 1;
simulink_opts.outputs.CPU_time = 0;
simulink_opts.outputs.sqp_iter = 1;
simulink_opts.outputs.x1 = 0;


%% create ocp solver
ocp = acados_ocp(ocp_model, ocp_opts, simulink_opts);

% solver initial guess
x_traj_init = rand(nx, N+1);
u_traj_init = zeros(nu, N);

%% call ocp solver
% update initial state
ocp.set('constr_x0', x0);

% set trajectory initialization
ocp.set('init_x', x_traj_init); % states
ocp.set('init_u', u_traj_init); % inputs
ocp.set('init_pi', zeros(nx, N)); % multipliers for dynamics equality constraints

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

%% simulink test
cd c_generated_code
make_sfun; % ocp solver

%% Run Simulink example block
out_sim = sim('initialization_test_simulink', 'SaveOutput', 'on');
disp('successfully ran simulink_model_advanced_closed_loop');

%% Evaluation
fprintf('\nTest results on SIMULINK simulation.\n')

disp('checking KKT residual')
kkt_signal = out_sim.logsout.getElement('KKT_residual');
if any(kkt_signal.Values.data > 1e-6)
    disp('failed');
    quit(1);
end

sqp_iter_signal = out_sim.logsout.getElement('sqp_iter');
disp('checking SQP iter, QP should take 1 SQP iter.')
if any(sqp_iter_signal.Values.Data ~= 1)
    disp('failed');
    quit(1);
end

status_signal = out_sim.logsout.getElement('status');
disp('checking status.')
if any(status_signal.Values.Data)
    disp('failed');
    quit(1);
end

utraj_signal = out_sim.logsout.getElement('utraj');
u_simulink = utraj_signal.Values.Data(1, :);
disp('checking u values.')
if any((u_simulink(:) - utraj(:)) > 1e-8)
    disp('failed');
    quit(1);
end

xtraj_signal = out_sim.logsout.getElement('xtraj');
x_simulink = xtraj_signal.Values.Data(1, :);
disp('checking x values.')
if any((x_simulink(:) - xtraj(:)) > 1e-8)
    disp('failed');
    quit(1);
end

