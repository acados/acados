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


clear all; clc;

check_acados_requirements()

% initial state
x0 = [0; 0; 0; 0];  % start at stable position

%% OCP DESCRIPTION
ocp = AcadosOcp();

%% IVP DESCRIPTION
sim = AcadosSim();

%% SOLVER OPTIONS

%% discretization
h = 0.01; % sampling time = length of first shooting interval
N = 20; % number of shooting intervals
% nonuniform discretization
shooting_nodes = [0.0 0.01, 0.05*(1:N-1)];
T = shooting_nodes(end);

ocp.solver_options.tf = T;
ocp.solver_options.N_horizon = N;
ocp.solver_options.shooting_nodes = shooting_nodes;
ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
% FULL_CONDENSING_HPIPM, PARTIAL_CONDENSING_HPIPM
% FULL_CONDENSING_QPOASES, PARTIAL_CONDENSING_OSQP
ocp.solver_options.qp_solver_cond_N = 5; % for partial condensing
ocp.solver_options.globalization = 'MERIT_BACKTRACKING'; % turns on globalization

% we add some model-plant mismatch by choosing different integration
% methods for model (within the OCP) and plant:

% integrator model
model_integrator_type = 'ERK';
model_sim_method_num_stages = 1;
model_sim_method_num_steps = 2;

ocp.solver_options.sim_method_num_stages = model_sim_method_num_stages;
ocp.solver_options.sim_method_num_steps = model_sim_method_num_steps;
ocp.solver_options.integrator_type = model_integrator_type;

% integrator plant
plant_integrator_type = 'IRK';
plant_sim_method_num_stages = 3;
plant_sim_method_num_steps = 3;

sim.solver_options.num_stages = plant_sim_method_num_stages;
sim.solver_options.num_steps = plant_sim_method_num_steps;
sim.solver_options.Tsim = h;
sim.solver_options.integrator_type = plant_integrator_type;

%% MODEL with mass as parameter
model = get_pendulum_on_cart_model(h, true);
ocp.model = model;
sim.model = model;

% dimensions
nx = model.x.rows();
nu = model.u.rows();

%% COST: nonlinear-least squares cost
ocp.cost.cost_type_0 = 'NONLINEAR_LS';
ocp.cost.cost_type = 'NONLINEAR_LS';
ocp.cost.cost_type_e = 'NONLINEAR_LS';

W_x = diag([1e3, 1e3, 1e-2, 1e-2]);
W_u = 1e-2;

model.cost_y_expr_0 = model.u;
model.cost_y_expr = vertcat(model.x, model.u);
model.cost_y_expr_e = model.x;

ocp.cost.W_0 = W_u;
ocp.cost.W = blkdiag(W_x, W_u);
ocp.cost.W_e = W_x;

% initialize reference to zero, can be changed after solver creation
ocp.cost.yref_0 = zeros(size(model.cost_y_expr_0));
ocp.cost.yref = zeros(size(model.cost_y_expr));
ocp.cost.yref_e = zeros(size(model.cost_y_expr_e));

%% CONSTRAINTS

U_max = 80;
ocp.constraints.constr_type = 'AUTO';
ocp.constraints.constr_type_0 = 'AUTO';

model.con_h_expr_0 = model.u;
ocp.constraints.lh_0 = -U_max;
ocp.constraints.uh_0 = U_max;

model.con_h_expr = model.u;
ocp.constraints.lh = -U_max;
ocp.constraints.uh = U_max;

ocp.constraints.x0 = x0;

%% OCP SOLVER
ocp_solver = AcadosOcpSolver(ocp);

% set parameter for all stages
for i = 0:N
    ocp_solver.set('p', 1., i);
end

%% SIM SOLVER/INTEGRATOR
sim_solver = AcadosSimSolver(sim);

% set parameter
sim_solver.set('p', 1.05); % model-plant mismatch in the parameters

%% SIMULATION
N_sim = 150;

% preallocate memory
x_sim = zeros(nx, N_sim+1);
u_sim = zeros(nu, N_sim);

x_sim(:,1) = x0;

% time-variant reference: move the cart with constant velocity while
% keeping the pendulum in upwards position
v_mean = 1;

yref_0 = zeros(nu, 1);
yref = zeros(nx+nu, 1);
yref_e = zeros(nx, 1);

yref(3) = v_mean;
yref_e(3) = v_mean;

for i=1:N_sim
    % update initial state
    x0 = x_sim(:,i);
    ocp_solver.set('constr_x0', x0);

    % compute reference position on the nonuniform grid
    t = (i-1)*h;
    p_ref = (t + shooting_nodes)*v_mean;

    for k=1:N-1 % intermediate stages
        yref(1) = p_ref(k);
        ocp_solver.set('cost_y_ref', yref, k); % last argument is the stage
    end
    yref_e(1) = p_ref(k+1); % terminal stage
    ocp_solver.set('cost_y_ref_e', yref_e, N);

    % solve
    ocp_solver.solve();

    % get solution
    u0 = ocp_solver.get('u', 0);
    status = ocp_solver.get('status'); % 0 - success

    % simulate one step
    x_sim(:,i+1) = sim_solver.simulate(x0, u0);
    u_sim(:,i) = u0;
end


%% plots
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
