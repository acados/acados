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

% initial state
xcurrent = [0.2; 0; 0; 0];

%% discretization
N = 20; % number of shooting intervals
% nonuniform discretization
T = 1.0;
shooting_nodes = linspace(0, T, N+1);
h = T/N; % sampling time = length of first shooting interval

nlp_solver = 'sqp'; % sqp, sqp_rti
qp_solver = 'partial_condensing_hpipm';
% full_condensing_hpipm, partial_condensing_hpipm, full_condensing_qpoases
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

% % nonlinear-least squares cost
% ocp_model.set('cost_type', 'nonlinear_ls');
% ocp_model.set('cost_type_e', 'nonlinear_ls');
% 
% ocp_model.set('cost_expr_y', model.cost_expr_y);
% ocp_model.set('cost_expr_y_e', model.cost_expr_y_e);
% 
% W_x = diag([1e2, 1e2, 1e-2, 1e-2]);
% W_u = 1e-3;
% W = blkdiag(W_x, W_u);
% ocp_model.set('cost_W', W);
% ocp_model.set('cost_W_e', model.W_e);

% % external cost -> with detection linear least squares
% cost
ocp_model.set('cost_expr_ext_cost', model.expr_ext_cost);
ocp_model.set('cost_expr_ext_cost_e', model.expr_ext_cost_e);


% dynamics
ocp_model.set('dyn_type', 'explicit');
ocp_model.set('dyn_expr_f', model.expr_f_expl);

% constraints
ocp_model.set('constr_type', 'auto');
ocp_model.set('constr_expr_h', model.expr_h);
U_max = 80;
ocp_model.set('constr_lh', -U_max); % lower bound on h
ocp_model.set('constr_uh', U_max);  % upper bound on h
ocp_model.set('constr_x0', xcurrent);

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
N_sim = 100;

x_sim = zeros(nx, N_sim+1);
u_sim = zeros(nu, N_sim);

x_sim(:,1) = xcurrent;

for i=1:N_sim
    % update initial state
    xcurrent = x_sim(:,i);
    ocp.set('constr_x0', xcurrent);

    if i == 1 || i == floor(N_sim/2)
        % solve
        ocp.solve();
        % get solution
        u0 = ocp.get('u', 0);
        status = ocp.get('status'); % 0 - success
        x_lin = xcurrent;
        u_lin = u0;
        
        sens_u = zeros(nx, N);
        % get sensitivities w.r.t. initial state value with index
        field = 'ex'; % equality constraint on states
        stage = 0;
        for index = 0:nx-1
            ocp.eval_param_sens(field, stage, index);
            sens_u(index+1,:) = ocp.get('sens_u');
        end
    else
        % use feedback policy
        delta_x = xcurrent-x_lin;
        u0 = u_lin + sens_u(:, 1)' * delta_x;
    end

    % set initial state
    sim.set('x', xcurrent);
    sim.set('u', u0);

    % solve
    sim.solve();

    % get simulated state
    x_sim(:,i+1) = sim.get('xn');
    u_sim(:,i) = u0;
end

disp('final state')
format long e
disp(x_sim(:,N_sim+1))
%      1.392073955008204e-03
%      4.247720422461933e-05
%     -7.518679517751918e-05
%     -1.671811407900214e-04

%% Plots
ts = linspace(0, N_sim*h, N_sim+1);
figure; hold on;
States = {'p', 'theta', 'v', 'dtheta'};
v_mean = 0;
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
stairs(ts(1:end), [u_sim'; u_sim(end)])
ylabel('F [N]')
xlabel('t [s]')
grid on
