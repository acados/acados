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
clear all

model_path = fullfile(pwd,'..','pendulum_on_cart_model');
addpath(model_path)

check_acados_requirements()

% initial state
x0 = [0; 0; 0; 0];  % start at stable position

%% discretization
h = 0.01; % sampling time = length of first shooting interval
N = 20; % number of shooting intervals
% nonuniform discretization
shooting_nodes = [0.0 0.01, 0.05*(1:N-1)];
T = shooting_nodes(end);

nlp_solver = 'sqp'; % sqp, sqp_rti
qp_solver = 'partial_condensing_hpipm';
% full_condensing_hpipm, partial_condensing_hpipm, full_condensing_qpoases, full_condensing_daqp
qp_solver_cond_N = 5; % for partial condensing

% we add some model-plant mismatch by choosing different integration
% methods for model (within the OCP) and plant:

% integrator model
model_sim_method = 'erk';
model_sim_method_num_stages = 1;
model_sim_method_num_steps = 2;

% integrator plant
plant_sim_method = 'irk';
plant_sim_method_num_stages = 3;
plant_sim_method_num_steps = 3;

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

% nonlinear-least squares cost
ocp_model.set('cost_type', 'nonlinear_ls');
ocp_model.set('cost_type_e', 'nonlinear_ls');

ocp_model.set('cost_expr_y', model.cost_expr_y);
ocp_model.set('cost_expr_y_e', model.cost_expr_y_e);

W_x = diag([1e2, 1e2, 1e-2, 1e-2]);
W_u = 1e-3;
W = blkdiag(W_x, W_u);
ocp_model.set('cost_W', W);
ocp_model.set('cost_W_e', model.W_e);

% dynamics
ocp_model.set('dyn_type', 'explicit');
ocp_model.set('dyn_expr_f', model.expr_f_expl);

% constraints
ocp_model.set('constr_type', 'auto');
ocp_model.set('constr_expr_h', model.expr_h);
U_max = 80;
ocp_model.set('constr_lh', -U_max); % lower bound on h
ocp_model.set('constr_uh', U_max);  % upper bound on h
ocp_model.set('constr_x0', x0);

%% acados ocp set opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('shooting_nodes', shooting_nodes);

ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('sim_method', model_sim_method);
ocp_opts.set('sim_method_num_stages', model_sim_method_num_stages);
ocp_opts.set('sim_method_num_steps', model_sim_method_num_steps);

ocp_opts.set('qp_solver', qp_solver);
ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
% ... see ocp_opts.opts_struct to see what other fields can be set

%% create ocp solver
ocp = acados_ocp(ocp_model, ocp_opts);

x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);


%% plant: create acados integrator
% acados sim model
sim_model = acados_sim_model();
sim_model.set('name', [model_name '_plant']);
sim_model.set('T', h);

sim_model.set('sym_x', model.sym_x);
sim_model.set('sym_u', model.sym_u);
sim_model.set('sym_xdot', model.sym_xdot);
sim_model.set('dyn_type', 'implicit');
sim_model.set('dyn_expr_f', model.expr_f_impl);

% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('method', plant_sim_method);
sim_opts.set('num_stages', plant_sim_method_num_stages);
sim_opts.set('num_steps', plant_sim_method_num_steps);

sim = acados_sim(sim_model, sim_opts);

%% Simulation
N_sim = 150;

x_sim = zeros(nx, N_sim+1);
u_sim = zeros(nu, N_sim);

x_sim(:,1) = x0;

% time-variant reference: move the cart with constant velocity while
% keeping the pendulum in upwards position
v_mean = 1;
yref = zeros(nx+nu, 1);
yref_e = zeros(nx, 1);

yref(3) = v_mean;
yref_e(3) = v_mean;

for i=1:N_sim
    % update initial state
    x0 = x_sim(:,i);
    ocp.set('constr_x0', x0);

    % compute reference position on the nonuniform grid
    t = (i-1)*h;
    p_ref = (t + shooting_nodes)*v_mean;

    for k=0:N-1
        yref(1) = p_ref(k+1);
        ocp.set('cost_y_ref', yref, k);
    end
    yref_e(1) = p_ref(k+1);
    ocp.set('cost_y_ref_e', yref_e, N);

    % solve
    ocp.solve();

    % get solution
    u0 = ocp.get('u', 0);
    status = ocp.get('status'); % 0 - success

    % set initial state
    sim.set('x', x0);
    sim.set('u', u0);

    % solve
    sim_status = sim.solve();
    if sim_status ~= 0
        disp(['acados integrator returned error status ', num2str(sim_status)])
    end

    % get simulated state
    x_sim(:,i+1) = sim.get('xn');
    u_sim(:,i) = u0;
end


%% Plots
ts = linspace(0, N_sim*h, N_sim+1);
figure; hold on;
States = {'p', 'theta', 'v', 'dtheta'};
p_ref = ts*v_mean;

y_ref = zeros(nx, N_sim+1);
y_ref(1, :) = p_ref;
y_ref(3, :) = v_mean;

for i=1:length(States)
    subplot(length(States), 1, i);
    grid on; hold on;
    plot(ts, x_sim(i,:)); 
    plot(ts, y_ref(i, :)); 
    ylabel(States{i});
    xlabel('t [s]')
    legend('closed-loop', 'reference')
end

figure
stairs(ts, [u_sim'; u_sim(end)])
ylabel('F [N]')
xlabel('t [s]')
grid on
