%
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

if 1
	compile_interface = 'auto';
	codgen_model = 'true';
	gnsf_detect_struct = 'true';
else
	compile_interface = 'auto';
	codgen_model = 'false';
	gnsf_detect_struct = 'false';
end

nlp_solver = 'sqp';
%nlp_solver = 'sqp_rti';
nlp_solver_exact_hessian = 'false';
%nlp_solver_exact_hessian = 'true';
regularize_method = 'no_regularize';
%regularize_method = 'project';
%regularize_method = 'project_reduc_hess';
%regularize_method = 'mirror';
%regularize_method = 'convexify';
nlp_solver_max_iter = 1000;
nlp_solver_tol_stat = 1e-6;
nlp_solver_tol_eq   = 1e-6;
nlp_solver_tol_ineq = 1e-6;
nlp_solver_tol_comp = 1e-6;
nlp_solver_step_length = 0.2;
nlp_solver_ext_qp_res = 1; % with 10 nothing changes
qp_solver = 'partial_condensing_hpipm';
%qp_solver = 'full_condensing_hpipm';
%qp_solver = 'full_condensing_qpoases';
qp_solver_cond_N = nb_steps/2; %5;
qp_solver_cond_ric_alg = 0;
qp_solver_ric_alg = 0;
qp_solver_warm_start = 0;
%sim_method = 'erk';
%sim_method = 'irk';
sim_method = 'irk_gnsf';
sim_method_num_stages = 4;
sim_method_num_steps = 3;
%cost_type = 'linear_ls';
cost_type = 'nonlinear_ls';
%cost_type = 'ext_cost';
model_name = 'ocp_swarming';


%% Model

model = swarming_model(S);

% Dimensions
nx = model.nx;
nu = model.nu;
ny = model.ny; % number of outputs in lagrange term
ny_e = model.ny_e; % number of outputs in mayer term

nbx = 0;
nbu = 0;
ng = 0;
ng_e = 0;
nh = nu;
nh_e = 0;

% Cost
W = eye(ny); % weight matrix in lagrange term
W_e = eye(ny_e); % weight matrix in mayer term

y_ref = zeros(ny, 1); % output reference in lagrange term
y_ref_e = zeros(ny_e,1); % output reference in mayer term

% Constraints
% x0 = [S.Pos0(:); S.Vel0(:)];
%rand('seed', 1);
pos0 = 10*rand(3*N,1);
vel0 = 2*rand(3*N,1);
x0 = [pos0; vel0];

lh = - max_a * ones(nh, 1);
uh = max_a * ones(nh, 1);
%lh_e = zeros(nh_e, 1);
%uh_e = zeros(nh_e, 1);

%% Acados ocp model
ocp_model = acados_ocp_model();
ocp_model.set('name', model_name);
ocp_model.set('T', T);

% Symbolics
ocp_model.set('sym_x', model.sym_x);
ocp_model.set('sym_u', model.sym_u);
ocp_model.set('sym_xdot', model.sym_xdot);

% Cost
ocp_model.set('cost_type', cost_type);
ocp_model.set('cost_type_e', cost_type);

if strcmp(cost_type, 'nonlinear_ls')
	ocp_model.set('cost_expr_y', model.expr_y);
	ocp_model.set('cost_expr_y_e', model.expr_y_e);
	ocp_model.set('cost_W', W);
	ocp_model.set('cost_W_e', W_e);
	ocp_model.set('cost_y_ref', y_ref);
	ocp_model.set('cost_y_ref_e', y_ref_e);
else % ext_cost
	ocp_model.set('cost_expr_ext_cost', model.expr_ext_cost);
	ocp_model.set('cost_expr_ext_cost_e', model.expr_ext_cost);
end

% Dynamics
if (strcmp(sim_method, 'erk'))
	ocp_model.set('dyn_type', 'explicit');
	ocp_model.set('dyn_expr_f', model.expr_f_expl);
elseif (strcmp(sim_method, 'irk') | strcmp(sim_method, 'irk_gnsf'))
	ocp_model.set('dyn_type', 'implicit');
	ocp_model.set('dyn_expr_f', model.expr_f_impl);
else
	ocp_model.set('dyn_type', 'discrete');
	ocp_model.set('dyn_expr_phi', model.expr_phi);
end

% Constraints
ocp_model.set('constr_x0', x0);

ocp_model.set('constr_expr_h', model.expr_h);
ocp_model.set('constr_lh', lh);
ocp_model.set('constr_uh', uh);
% ocp_model.set('constr_expr_h_e', model.expr_h_e);
% ocp_model.set('constr_lh_e', lh_e);
% ocp_model.set('constr_uh_e', uh_e);

ocp_model.model_struct

%% Acados ocp options

ocp_opts = acados_ocp_opts();
ocp_opts.set('compile_interface', compile_interface);
ocp_opts.set('codgen_model', codgen_model);
ocp_opts.set('param_scheme_N', nb_steps);
if (exist('shooting_nodes', 'var'))
	ocp_opts.set('shooting_nodes', shooting_nodes);
end
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('nlp_solver_exact_hessian', nlp_solver_exact_hessian);
ocp_opts.set('regularize_method', regularize_method);
ocp_opts.set('nlp_solver_ext_qp_res', nlp_solver_ext_qp_res);
if (strcmp(nlp_solver, 'sqp'))
	ocp_opts.set('nlp_solver_max_iter', nlp_solver_max_iter);
	ocp_opts.set('nlp_solver_tol_stat', nlp_solver_tol_stat);
	ocp_opts.set('nlp_solver_tol_eq', nlp_solver_tol_eq);
	ocp_opts.set('nlp_solver_tol_ineq', nlp_solver_tol_ineq);
	ocp_opts.set('nlp_solver_tol_comp', nlp_solver_tol_comp);
    ocp_opts.set('nlp_solver_step_length', nlp_solver_step_length);
    ocp_opts.set('qp_solver_tol_stat', .1 * nlp_solver_tol_stat);
end
ocp_opts.set('qp_solver', qp_solver);
if (strcmp(qp_solver, 'partial_condensing_hpipm'))
	ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
	ocp_opts.set('qp_solver_ric_alg', qp_solver_ric_alg);
end
ocp_opts.set('qp_solver_cond_ric_alg', qp_solver_cond_ric_alg);
ocp_opts.set('qp_solver_warm_start', qp_solver_warm_start);
ocp_opts.set('sim_method', sim_method);
ocp_opts.set('sim_method_num_stages', sim_method_num_stages);
ocp_opts.set('sim_method_num_steps', sim_method_num_steps);
if (strcmp(sim_method, 'irk_gnsf'))
	ocp_opts.set('gnsf_detect_struct', gnsf_detect_struct);
end

ocp_opts.opts_struct

%% Acados ocp

% Create ocp
ocp = acados_ocp(ocp_model, ocp_opts);
ocp
ocp.C_ocp

% Set trajectory initialization
step_mat = repmat((0:1:nb_steps),3*N,1);
pos0_traj = repmat(pos0,1,nb_steps+1) + v_ref*dt*repmat(u_ref,N,nb_steps+1).*step_mat;
x_traj_init = [pos0_traj; ...
    v_ref*repmat(u_ref,N,nb_steps+1)];
u_traj_init = zeros(nu, nb_steps);
ocp.set('init_x', x_traj_init);
ocp.set('init_u', u_traj_init);

% Solve
tic;
ocp.solve();
simulation_time = toc;
disp(strcat('Simulation time: ',num2str(simulation_time)));

%x0(1) = 1.5;
%ocp.set('constr_x0', x0);
%ocp.set('cost_y_ref', 1);
% If not set, the trajectory is initialized with the previous solution

% Get solution
u = ocp.get('u');
x = ocp.get('x');

%% Statistics

status = ocp.get('status');
sqp_iter = ocp.get('sqp_iter');
time_tot = ocp.get('time_tot');
time_lin = ocp.get('time_lin');
time_reg = ocp.get('time_reg');
time_qp_sol = ocp.get('time_qp_sol');

fprintf('\nstatus = %d, sqp_iter = %d,  time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n', status, sqp_iter, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3);
% time_ext = %f [ms], time_ext*1e3,

ocp.print('stat');

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

if (strcmp(nlp_solver, 'sqp'))
	figure;
    stat = ocp.get('stat');
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
