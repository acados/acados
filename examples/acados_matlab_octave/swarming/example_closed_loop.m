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

clear all

% Check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	disp('ERROR: env.sh has not been sourced! Before executing this example, run:');
	disp('source env.sh');
	return;
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
rand('seed', 2);
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


%% Acados sim model
sim = AcadosSim();
sim.model = model;
sim.solver_options.Tsim = dt;
sim.solver_options.integrator_type = 'ERK';  % 'ERK', 'IRK'
sim.solver_options.num_stages = 4;
sim.solver_options.num_steps = 3;
sim.solver_options.sens_forw = true;
sim.solver_options.compile_interface = 'AUTO';

%% Acados simulation

% Create sim
sim_solver = AcadosSimSolver(sim);


%% Closed loop simulation

T_sim = 20; % time of the whole simulation
nb_steps_sim = floor(T_sim/dt); % time of time steps during the simulation
x_history = zeros(nx, nb_steps_sim+1);
x_history(:,1) = x0;
u_history = zeros(nu, nb_steps_sim);

% Set state and input trajectory initialization
step_mat = repmat((0:1:nb_steps),3*N,1);
pos0_traj = repmat(pos0,1,nb_steps+1) + v_ref*dt*repmat(u_ref,N,nb_steps+1).*step_mat;
x_traj_init = [pos0_traj; ...
    v_ref*repmat(u_ref,N,nb_steps+1)];
u_traj_init = zeros(nu, nb_steps);

% Initialize variables
status = zeros(1, nb_steps_sim);
sqp_iter = zeros(1, nb_steps_sim);
time_tot = zeros(1, nb_steps_sim);
time_lin = zeros(1, nb_steps_sim);
time_qp_sol = zeros(1, nb_steps_sim);

tic;

for k = 1:nb_steps_sim

	% Set initial condition x0
	ocp_solver.set('constr_x0', x_history(:,k));
%     ocp_solver.set('constr_expr_h', model.expr_h);
%     ocp_solver.set('constr_lh', lh);
%     ocp_solver.set('constr_uh', uh);

	% Set trajectory initialization (if not, set internally using previous solution)
	ocp_solver.set('init_x', x_traj_init);
	ocp_solver.set('init_u', u_traj_init);

	% solve OCP
	ocp_solver.solve();

    status(k) = ocp_solver.get('status');
    sqp_iter(k) = ocp_solver.get('sqp_iter');
    time_tot(k) = ocp_solver.get('time_tot');
    time_lin(k) = ocp_solver.get('time_lin');
    time_qp_sol(k) = ocp_solver.get('time_qp_sol');

    fprintf('\nstatus = %d, sqp_iter = %d, time_tot = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms])\n', status(k), sqp_iter(k), time_tot(k)*1e3, time_lin(k)*1e3, time_qp_sol(k)*1e3);

	% Get solution for initialization of next NLP
	x_traj = ocp_solver.get('x');
	u_traj = ocp_solver.get('u');

	% Shift trajectory for initialization
	x_traj_init = [x_traj(:,2:end), x_traj(:,end)];

	% Get solution for simulation
	u_history(:,k) = ocp_solver.get('u', 0);

	% Set initial state of simulation
	sim_solver.set('x', x_history(:,k));
	% set input in sim
	sim_solver.set('u', u_history(:,k));

	% Simulate state
	sim_solver.solve();

	% Get new state
	x_history(:,k+1) = sim_solver.get('xn');

end

simulation_time = toc;
disp(strcat('Simulation time: ',num2str(simulation_time)));

%% Extract trajectories

fontsize = 12;

time_history = linspace(0,T_sim,nb_steps_sim+1)';
x_history = x_history';
u_history = u_history';
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
plot(time_history(1:(end-1)), u_history);
xlim([0 time_history(end-1)]);
xlabel('Time [s]','fontsize',fontsize);
ylabel('Control inputs [m/s^2]','fontsize',fontsize);

%% Show solver convergence

figure;
plot([1: nb_steps_sim], (sqp_iter), 'r-x');
xlabel('Iteration','fontsize',fontsize)
ylabel('Nb SQP iteration','fontsize',fontsize);
ylim([0 Inf]);

figure;
plot([1: nb_steps_sim], (time_tot*1e3), 'b-x');
xlabel('Iteration','fontsize',fontsize)
ylabel('Simulation time [ms]','fontsize',fontsize);
ylim([0 Inf]);


if any(status)== 0
	fprintf('\nsuccess!\n\n');
else
	fprintf('\nsolution failed in some steps!\n\n');
end

if is_octave
    waitforbuttonpress;
end