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

%% Test of native matlab interface
clear all; close all; clc;
check_acados_requirements()

% discretization
h = 0.01;
N_horizon = 100;
T = N_horizon*h;

% initial state
x0 = [0.1; 0; 0; 0];

model = F8_crusader_model();

% dimension
nx = length(model.x);
nu = length(model.u);
ny = nx + nu;                  % number of outputs in lagrange term
ny_e = nx;                     % number of outputs in mayer term
nbx = nx;                      % number of state bounds
nbu = nu;                      % number of input bounds

% setup OCP
ocp = AcadosOcp();
ocp.model = model;


% integrator
ocp.solver_options.integrator_type = 'ERK';
ocp.solver_options.sim_method_num_stages = 4;
ocp.solver_options.sim_method_num_steps = 1;

ocp.solver_options.tf = T;
ocp.solver_options.N_horizon = N_horizon;

% NLP solver
ocp.solver_options.nlp_solver_type = 'SQP_RTI';

% QP solver
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_iter_max = 50;
ocp.solver_options.qp_solver_cond_N = 5; % number of shooting nodes after partial_condensing
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';

% cost
ocp.cost.cost_type = 'LINEAR_LS';
ocp.cost.cost_type_e = 'LINEAR_LS';

Vx = zeros(ny,nx); Vx(1:nx,:) = eye(nx);        % state-to-output matrix in lagrange term
Vu = zeros(ny,nu); Vu(nx+1:ny,:) = eye(nu);     % input-to-output matrix in lagrange term
Vx_e = zeros(ny_e,nx); Vx_e(1:nx,:) = eye(nx);  % state-to-output matrix in mayer term
W = diag([10, 0.1, 0, 1, 0.01]);                % cost weights in lagrange term
W_e = W(1:ny_e,1:ny_e);                         % cost weights in mayer term
yref = zeros(ny,1);                             % references
yref_e = zeros(ny_e,1);

ocp.cost.Vx = Vx;
ocp.cost.Vu = Vu;
ocp.cost.Vx_e = Vx_e;

ocp.cost.W = W;
ocp.cost.W_e = W_e;

ocp.cost.yref = yref;
ocp.cost.yref_e = yref_e;

% constraints
ocp.constraints.x0 = x0;

ocp.constraints.idxbu = (1:nu) - 1;
ocp.constraints.lbu = -0.3;
ocp.constraints.ubu = 0.5;

ocp.constraints.idxbx = (1:nx) - 1;
ocp.constraints.lbx = [-0.2; -1; -1; -0.3];
ocp.constraints.ubx = [0.4; 1; 1; 0.5];

% create OCP solver
ocp_solver = AcadosOcpSolver(ocp);

%% Simulation
simulation_time = 10;  % [s]
N_sim = simulation_time/h;

% initialize data structs
x_sim = zeros(nx, N_sim+1);
u_sim = zeros(nu, N_sim+1);
cost_sim = zeros(1, N_sim+1);

x_sim(:, 1) = x0;
u_sim(:, 1) = zeros(nu, 1);
cost_sim(1, 1) = 0;

% set trajectory initialization
ocp_solver.set('init_x', x0 * ones(1,N_horizon+1));
ocp_solver.set('init_u', zeros(nu, N_horizon));

% time-varying reference trajectory
x1ref_FUN = @(t) 0.4.*(-(0.5./(1+exp(t./0.1-0.8))) + (1./(1+exp(t./0.1-30))) - 0.4);
t = 0:h:simulation_time;
x1ref = 0.4.*(-(0.5./(1+exp(t./0.1-0.8))) + (1./(1+exp(t./0.1-30))) - 0.4);

% run mpc
fprintf('Simulation started.  It might take a while...\n')
tic;
for i = 1:N_sim

    % update reference (full preview)
    t_ref = (i-1:i+N_horizon).*h;
    x1_ref = x1ref_FUN(t_ref);
    for j = 0:N_horizon-1
        yref(1) = x1_ref(j+1);
        ocp_solver.set('cost_y_ref', yref, j);
    end
    yref_e(1) = x1_ref(N_horizon+1);
    ocp_solver.set('cost_y_ref_e', yref_e, N_horizon);

    % solve ocp
    ocp_solver.solve();
    status = ocp_solver.get('status');      % 0 - success
    if status ~= 0
        error('acados returned status %d in closed loop iteration %d. Exiting.', status, i);
    end

    % get solution at initial shooting node
    x0 = ocp_solver.get('x', 0);
    u0 = ocp_solver.get('u', 0);
    x_sim(:, i+1) = x0;
    u_sim(:, i+1) = u0;
    cost_sim(1, i+1) = ocp_solver.get_cost();

    % update initial state
    % TODO use an integrator here instead of the OCP state (we use SQP_RTI)
    x0 = ocp_solver.get('x', 1);
    ocp_solver.set('constr_x0', x0);

end
tElapsed = toc
fprintf('Simulation finished!\n')

%% Plot
figure; hold on; grid on;
plot(t, x_sim, t, x1ref, '--')
legend('x1', 'x2', 'x3', 'u', 'x1Ref')
ylabel('(augmented) state')
xlabel('time [s]')

figure; hold on; grid on;
plot(t, u_sim)
plot(t, ocp.constraints.lbu, 'k--')
plot(t, ocp.constraints.ubu, 'k--')
legend({'udot'})
ylabel('control input rate')
xlabel('time [s]')

figure; hold on; grid on;
plot(t, cost_sim)
legend('the cost curve')
ylabel('cost')
xlabel('time [s]')
