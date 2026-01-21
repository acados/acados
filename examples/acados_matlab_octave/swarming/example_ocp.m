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

% Author: Enrica

% This file allows the control of a swarm of robots. The swarm is composed
% by N agents with decoupled, linear dynamics. The goal is to achieve
% coordinated motion from random position and velocities.
%


%% Test of native matlab interface

clear all;
close all;

% Check that env.sh has been runz2
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end

%% Arguments

% Time parameters
dt = 0.1; % discretization step
T = 8; % total time horizon of the simulation
nb_steps = floor(T/dt); % nb of time steps along the simulation

% Structure S with the swarming parameters
S.N = 3; % number of agents in the swarm
S.d_ref = 5; % reference distance among every couple of neighboring agents
S.u_ref = [1;0;0]; % reference direction of velocity for all agents
S.v_ref = 6; % reference speed for all agents
S.max_a = 2;

% Rename swarming parameters
N = S.N;
u_ref = S.u_ref;
v_ref = S.v_ref;
max_a = S.max_a;

%% Model
model = get_swarming_model(S);
% Dimensions
nx = length(model.x);
nu = length(model.u);

nbx = 0;
nbu = 0;
ng = 0;
ng_e = 0;
nh = nu;
nh_e = 0;
%% Acados ocp
ocp = AcadosOcp();
ocp.model = model;

tol = 1e-6;
ocp.solver_options.tf = T;
ocp.solver_options.N_horizon = nb_steps;
ocp.solver_options.nlp_solver_type = 'SQP'; % 'SQP_RTI'
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'; % 'EXACT', 'GAUSS_NEWTON'
ocp.solver_options.regularize_method = 'NO_REGULARIZE';% NO_REGULARIZE, PROJECT, PROOJECT_REDUC_HESS, MIRROR, CONVEXIFY
ocp.solver_options.nlp_solver_max_iter = 1000;
ocp.solver_options.nlp_solver_tol_stat = tol;
ocp.solver_options.nlp_solver_tol_eq = tol;
ocp.solver_options.nlp_solver_tol_ineq = tol;
ocp.solver_options.nlp_solver_tol_comp = tol;
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_cond_N = nb_steps/2; % for partial condensing
ocp.solver_options.nlp_solver_step_length = 0.2;
ocp.solver_options.qp_solver_cond_ric_alg = 0;
ocp.solver_options.qp_solver_ric_alg = 0;
ocp.solver_options.qp_solver_warm_start = 0;
ocp.solver_options.integrator_type = 'ERK';
ocp.solver_options.nlp_solver_ext_qp_res = 1;
ocp.solver_options.sim_method_num_stages = 4;
ocp.solver_options.sim_method_num_steps = 3;
ocp.solver_options.compile_interface = 'AUTO';
ocp.model.name ='ocp_swarming';

% Cost (the cost_y_expr is designed in get_swarming_model.py)
ocp.cost.cost_type = 'NONLINEAR_LS';
ocp.cost.cost_type_e = 'NONLINEAR_LS';

ny = length(ocp.model.cost_y_expr);
ny_e = length(ocp.model.cost_y_expr_e);

W = eye(ny); % weight matrix in lagrange term
W_e = eye(ny_e); % weight matrix in mayer term

y_ref = zeros(ny, 1); % output reference in lagrange term
y_ref_e = zeros(ny_e,1); % output reference in mayer term

ocp.cost.W = W;
ocp.cost.W_e = W_e;
ocp.cost.yref = y_ref; 
ocp.cost.yref_e = y_ref_e;

% Constraints
expr_h = ocp.model.u; % constraints only on control inputs, for now
% rand('seed', 2);
pos0 = 10*rand(3*N,1);
vel0 = 2*rand(3*N,1);
x0 = [pos0; vel0];
lh = - max_a * ones(nh, 1);
uh = max_a * ones(nh, 1);
%lh_e = zeros(nh_e, 1);
%uh_e = zeros(nh_e, 1);
% expr_h_e = sym_x;

ocp.constraints.x0 = x0;
ocp.model.con_h_expr_0 = expr_h;
ocp.constraints.lh_0 = lh;
ocp.constraints.uh_0 = uh;
ocp.model.con_h_expr = expr_h;
ocp.constraints.lh = lh;
ocp.constraints.uh = uh;
% ocp.model.con_h_expr_e = expr_h_e;
% ocp.constraints.lh_e = lh_e;
% ocp.constraints.uh_e = uh_e;

%% Acados ocp solver
% Create ocp
ocp_solver = AcadosOcpSolver(ocp);

% Set trajectory initialization
step_mat = repmat((0:1:nb_steps),3*N,1);
pos0_traj = repmat(pos0,1,nb_steps+1) + v_ref*dt*repmat(u_ref,N,nb_steps+1).*step_mat;
x_traj_init = [pos0_traj; ...
    v_ref*repmat(u_ref,N,nb_steps+1)];
u_traj_init = zeros(nu, nb_steps);
ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);

% Solve
tic;
ocp_solver.solve();
simulation_time = toc;
disp(strcat('Simulation time: ',num2str(simulation_time)));

%x0(1) = 1.5;
%ocp_solver.set('constr_x0', x0);
%ocp_solver.set('cost_y_ref', 1);
% If not set, the trajectory is initialized with the previous solution

% Get solution
u = ocp_solver.get('u');
x = ocp_solver.get('x');

%% Statistics

status = ocp_solver.get('status');
sqp_iter = ocp_solver.get('sqp_iter');
time_tot = ocp_solver.get('time_tot');
time_lin = ocp_solver.get('time_lin');
time_reg = ocp_solver.get('time_reg');
time_qp_sol = ocp_solver.get('time_qp_sol');

fprintf('\nstatus = %d, sqp_iter = %d,  time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n', status, sqp_iter, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3);
% time_ext = %f [ms], time_ext*1e3,

ocp_solver.print('stat');

%% Extract trajectories

fontsize = 12;

time_history = linspace(0,T,nb_steps+1)';
x_history = x';
u_history = u';
pos_history = x_history(:,1:3*N);
vel_history = x_history(:,(3*N+1):end);


%% Plots

% Plot trajectories of the agents
figure;
for agent = 1:N
    hold on;
    plot3(pos_history(:,(agent-1)*3+1), pos_history(:,(agent-1)*3+2), ...
        - pos_history(:,(agent-1)*3+3));
end
% title('Agents trajectories');
xlabel('X Position [m]','fontsize',fontsize);
ylabel('Y Position [m]','fontsize',fontsize);
zlabel('Z Position [m]','fontsize',fontsize);
view(2);

% Plot control inputs of the agents
figure;
plot(time_history(1:(end-1)), u);
xlim([0 time_history(end-1)]);
xlabel('Time [s]','fontsize',fontsize);
ylabel('Control inputs [m/s^2]','fontsize',fontsize);

%% Show solver convergence

if (strcmp(ocp.solver_options.nlp_solver_type, 'SQP'))
	figure;
    stat = ocp_solver.get('stat');
	plot([0: sqp_iter], log10(stat(:,2)), 'r-x');
	hold on
	plot([0: sqp_iter], log10(stat(:,3)), 'b-x');
	plot([0: sqp_iter], log10(stat(:,4)), 'g-x');
	plot([0: sqp_iter], log10(stat(:,5)), 'k-x');
	hold off
	xlabel('Iteration','fontsize',fontsize)
	ylabel('Residuals (log10)','fontsize',fontsize)
    legend('res g','res b','res d','res m','fontsize',fontsize)
end

if status == 0
	fprintf('\nsuccess!\n\n');
else
	fprintf('\nsolution failed!\n\n');
end

if is_octave
    waitforbuttonpress;
end
