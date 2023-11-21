%% Test of native matlab interface
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
dt = 0.1; % sample time
T = 4; % prediction horizon [s]
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

% Simulation
sim_method = 'erk'; % erk, irk, irk_gnsf
sim_sens_forw = 'false'; % true, false
sim_num_stages = 4;
sim_num_steps = 3;

% OCP
nlp_solver = 'sqp'; % sqp, sqp_rti
nlp_solver_exact_hessian = 'false';
regularize_method = 'no_regularize'; % no_regularize, project,...
	% project_reduc_hess, mirror, convexify
nlp_solver_max_iter = 1000;
nlp_solver_tol_stat = 1e-6;
nlp_solver_tol_eq   = 1e-6;
nlp_solver_tol_ineq = 1e-6;
nlp_solver_tol_comp = 1e-6;
nlp_solver_step_length = 0.2;
nlp_solver_ext_qp_res = 1; % with 10 nothing changes
qp_solver = 'partial_condensing_hpipm';
        % full_condensing_hpipm, partial_condensing_hpipm
qp_solver_cond_N = nb_steps/2;
qp_solver_warm_start = 0;
qp_solver_cond_ric_alg = 0; % 0: dont factorize hessian in the condensing; 1: factorize
qp_solver_ric_alg = 0; % HPIPM specific
ocp_sim_method = 'irk_gnsf'; % erk, irk, irk_gnsf
ocp_sim_method_num_stages = 4;
ocp_sim_method_num_steps = 3;
cost_type = 'nonlinear_ls'; % linear_ls, ext_cost

model_name = 'cl_swarming';

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
if (strcmp(ocp_sim_method, 'erk'))
	ocp_model.set('dyn_type', 'explicit');
	ocp_model.set('dyn_expr_f', model.expr_f_expl);
elseif (strcmp(ocp_sim_method, 'irk') | strcmp(ocp_sim_method, 'irk_gnsf'))
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
	ocp_opts.set('param_scheme_shooting_nodes', shooting_nodes);
end
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('nlp_solver_exact_hessian', nlp_solver_exact_hessian);
ocp_opts.set('regularize_method', regularize_method);
if (strcmp(nlp_solver, 'sqp'))
	ocp_opts.set('nlp_solver_max_iter', nlp_solver_max_iter);
	ocp_opts.set('nlp_solver_tol_stat', nlp_solver_tol_stat);
	ocp_opts.set('nlp_solver_tol_eq', nlp_solver_tol_eq);
	ocp_opts.set('nlp_solver_tol_ineq', nlp_solver_tol_ineq);
	ocp_opts.set('nlp_solver_tol_comp', nlp_solver_tol_comp);
    ocp_opts.set('nlp_solver_step_length', nlp_solver_step_length);
end
ocp_opts.set('qp_solver', qp_solver);
if (strcmp(qp_solver, 'partial_condensing_hpipm'))
	ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
	ocp_opts.set('qp_solver_ric_alg', qp_solver_ric_alg);
end
ocp_opts.set('qp_solver_cond_ric_alg', qp_solver_cond_ric_alg);
ocp_opts.set('qp_solver_warm_start', qp_solver_warm_start);
ocp_opts.set('sim_method', ocp_sim_method);
ocp_opts.set('sim_method_num_stages', ocp_sim_method_num_stages);
ocp_opts.set('sim_method_num_steps', ocp_sim_method_num_steps);
if (strcmp(ocp_sim_method, 'irk_gnsf'))
	ocp_opts.set('gnsf_detect_struct', gnsf_detect_struct);
end

ocp_opts.opts_struct

%% Acados ocp

% Create ocp
ocp = acados_ocp(ocp_model, ocp_opts);
ocp
ocp.C_ocp

%% Acados simulation model

sim_model = acados_sim_model();

% Symbolics
sim_model.set('sym_x', model.sym_x);
if isfield(model, 'sym_u')
	sim_model.set('sym_u', model.sym_u);
end
if isfield(model, 'sym_xdot')
	sim_model.set('sym_xdot', model.sym_xdot);
end
% model
sim_model.set('T', dt);
if (strcmp(sim_method, 'erk'))
	sim_model.set('dyn_type', 'explicit');
	sim_model.set('dyn_expr_f', model.expr_f_expl);
else % irk
	sim_model.set('dyn_type', 'implicit');
	sim_model.set('dyn_expr_f', model.expr_f_impl);
end

%sim_model.model_struct



%% Acados simulation options
sim_opts = acados_sim_opts();
sim_opts.set('compile_interface', compile_interface);
sim_opts.set('codgen_model', codgen_model);
sim_opts.set('num_stages', sim_num_stages);
sim_opts.set('num_steps', sim_num_steps);
sim_opts.set('method', sim_method);
sim_opts.set('sens_forw', sim_sens_forw);
if (strcmp(sim_method, 'irk_gnsf'))
	sim_opts.set('gnsf_detect_struct', gnsf_detect_struct);
end

%sim_opts.opts_struct



%% Acados simulation

% Create sim
sim = acados_sim(sim_model, sim_opts);


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
	ocp.set('constr_x0', x_history(:,k));
%     ocp.set('constr_expr_h', model.expr_h);
%     ocp.set('constr_lh', lh);
%     ocp.set('constr_uh', uh);

	% Set trajectory initialization (if not, set internally using previous solution)
	ocp.set('init_x', x_traj_init);
	ocp.set('init_u', u_traj_init);

	% solve OCP
	ocp.solve();

    status(k) = ocp.get('status');
    sqp_iter(k) = ocp.get('sqp_iter');
    time_tot(k) = ocp.get('time_tot');
    time_lin(k) = ocp.get('time_lin');
    time_qp_sol(k) = ocp.get('time_qp_sol');

    fprintf('\nstatus = %d, sqp_iter = %d, time_tot = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms])\n', status(k), sqp_iter(k), time_tot(k)*1e3, time_lin(k)*1e3, time_qp_sol(k)*1e3);

	% Get solution for initialization of next NLP
	x_traj = ocp.get('x');
	u_traj = ocp.get('u');
    
	% Shift trajectory for initialization
	x_traj_init = [x_traj(:,2:end), x_traj(:,end)];
	u_traj_init = [u_traj(:,2:end), u_traj(:,end)];

	% Get solution for simulation
	u_history(:,k) = ocp.get('u', 0);

	% Set initial state of simulation
	sim.set('x', x_history(:,k));
	% set input in sim
	sim.set('u', u_history(:,k));

	% Simulate state
	sim.solve();

	% Get new state
	x_history(:,k+1) = sim.get('xn');

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


if status == 0
	fprintf('\nsuccess!\n\n');
else
	fprintf('\nsolution failed!\n\n');
end

waitforbuttonpress;
return;
