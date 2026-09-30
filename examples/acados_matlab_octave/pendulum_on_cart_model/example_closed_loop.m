%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

clear all

model_name = 'ocp_pendulum';

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end


%% options
compile_interface = 'auto'; % true, false
% simulation
sim_method = 'irk'; % erk, irk, irk_gnsf
sim_sens_forw = 'false'; % true, false
sim_num_stages = 4;
sim_num_steps = 4;
% ocp
ocp_N = 100;
%nlp_solver = 'sqp_rti';
%nlp_solver_exact_hessian = 'false';
nlp_solver = 'sqp'; % sqp, sqp_rti
nlp_solver_exact_hessian = 'false';
regularize_method = 'project_reduc_hess'; % no_regularize, project,...
	% project_reduc_hess, mirror, convexify
%regularize_method = 'mirror';
%regularize_method = 'convexify';
nlp_solver_max_iter = 100;
qp_solver = 'partial_condensing_hpipm';
        % full_condensing_hpipm, partial_condensing_hpipm, full_condensing_qpoases
qp_solver_iter_max = 100;
qp_solver_cond_N = 5;
qp_solver_warm_start = 0;
qp_solver_cond_ric_alg = 0; % 0: dont factorize hessian in the condensing; 1: factorize
qp_solver_ric_alg = 0; % HPIPM specific
ocp_sim_method = 'erk'; % erk, irk, irk_gnsf
% ocp_sim_method = 'irk';
ocp_sim_method_num_stages = 4;
ocp_sim_method_num_steps = 1;
cost_type = 'linear_ls'; % linear_ls, ext_cost
% cost_type = 'ext_cost'; % linear_ls, ext_cost


%% create model entries
model = get_pendulum_on_cart_model();

h = 0.01;
T = ocp_N*h; % horizon length time

% dims
nx = length(model.x);
nu = length(model.u);

ny = nu+nx; % number of outputs in lagrange term
ny_e = nx; % number of outputs in mayer term

ng = 0; % number of general linear constraints intermediate stages
ng_e = 0; % number of general linear constraints final stage
nbx = 0; % number of bounds on state x


linear_constraints = 1; % 1: encode control bounds as bounds (efficient)
    % 0: encode control bounds as external CasADi functions
if linear_constraints
	nbu = nu;
	nh = 0;
	nh_e = 0;
else
	nbu = 0;
	nh = nu;
	nh_e = 0;
end

% cost
% linear least square cost: y^T * W * y, where y = Vx * x + Vu * u - y_ref
Vu = zeros(ny, nu); for ii=1:nu Vu(ii,ii)=1.0; end % input-to-output matrix in lagrange term
Vx = zeros(ny, nx); for ii=1:nx Vx(nu+ii,ii)=1.0; end % state-to-output matrix in lagrange term
Vx_e = zeros(ny_e, nx); for ii=1:nx Vx_e(ii,ii)=1.0; end % state-to-output matrix in mayer term
W = eye(ny); % weight matrix in lagrange term
for ii=1:nu W(ii,ii)=1e-2; end
for ii=nu+1:nu+nx/2 W(ii,ii)=1e3; end
for ii=nu+nx/2+1:nu+nx W(ii,ii)=1e-2; end
W_e = W(nu+1:nu+nx, nu+1:nu+nx); % weight matrix in mayer term
yr = zeros(ny, 1); % output reference in lagrange term
yr_e = zeros(ny_e, 1); % output reference in mayer term

% constraints
x0 = [0; pi; 0; 0];
%Jbx = zeros(nbx, nx); for ii=1:nbx Jbx(ii,ii)=1.0; end
%lbx = -4*ones(nbx, 1);
%ubx =  4*ones(nbx, 1);
Jbu = zeros(nbu, nu); for ii=1:nbu Jbu(ii,ii)=1.0; end
lbu = -80*ones(nu, 1);
ubu =  80*ones(nu, 1);



%% OCP and simulation formulations
ocp = AcadosOcp();
ocp.model = model;

if strcmp(cost_type, 'ext_cost')
	ocp.cost.cost_type_0 = 'EXTERNAL';
	ocp.cost.cost_type = 'EXTERNAL';
	ocp.cost.cost_type_e = 'EXTERNAL';
else
	ocp.cost.cost_type_0 = 'LINEAR_LS';
	ocp.cost.cost_type = 'LINEAR_LS';
	ocp.cost.cost_type_e = 'LINEAR_LS';
end
if strcmp(cost_type, 'linear_ls')
	ocp.cost.Vu_0 = Vu;
	ocp.cost.Vx_0 = Vx;
	ocp.cost.W_0 = W;
	ocp.cost.yref_0 = yr;
	ocp.cost.Vu = Vu;
	ocp.cost.Vx = Vx;
	ocp.cost.Vx_e = Vx_e;
	ocp.cost.W = W;
	ocp.cost.W_e = W_e;
	ocp.cost.yref = yr;
	ocp.cost.yref_e = yr_e;
else
	ocp.model.cost_expr_ext_cost_0 = 0.5 * model.u' * 1e-2 * model.u;
	ocp.model.cost_expr_ext_cost = 0.5 * model.x' * diag([1e3, 1e3, 1e-2, 1e-2]) * model.x + 0.5 * model.u' * 1e-2 * model.u;
	ocp.model.cost_expr_ext_cost_e = 0.5 * model.x' * diag([1e3, 1e3, 1e-2, 1e-2]) * model.x;
end

if strcmp(ocp_sim_method, 'erk')
	ocp.solver_options.integrator_type = 'ERK';
else
	ocp.solver_options.integrator_type = 'IRK';
end
ocp.constraints.x0 = x0;
if nh > 0
	ocp.model.con_h_expr_0 = model.u;
	ocp.constraints.lh_0 = lbu;
	ocp.constraints.uh_0 = ubu;
	ocp.model.con_h_expr = model.u;
	ocp.constraints.lh = lbu;
	ocp.constraints.uh = ubu;
else
	ocp.constraints.idxbu = (0:nu-1)';
	ocp.constraints.lbu = lbu;
	ocp.constraints.ubu = ubu;
end

ocp.solver_options.N_horizon = ocp_N;
ocp.solver_options.tf = T;
ocp.solver_options.nlp_solver_type = upper(nlp_solver);
if strcmp(nlp_solver_exact_hessian, 'true')
	ocp.solver_options.hessian_approx = 'EXACT';
else
	ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
end
ocp.solver_options.regularize_method = upper(regularize_method);
ocp.solver_options.nlp_solver_max_iter = nlp_solver_max_iter;
ocp.solver_options.qp_solver = upper(qp_solver);
ocp.solver_options.qp_solver_iter_max = qp_solver_iter_max;
if strcmp(qp_solver, 'partial_condensing_hpipm')
	ocp.solver_options.qp_solver_cond_N = qp_solver_cond_N;
	ocp.solver_options.qp_solver_cond_ric_alg = qp_solver_cond_ric_alg;
	ocp.solver_options.qp_solver_ric_alg = qp_solver_ric_alg;
	ocp.solver_options.qp_solver_warm_start = qp_solver_warm_start;
end
ocp.solver_options.sim_method_num_stages = ocp_sim_method_num_stages;
ocp.solver_options.sim_method_num_steps = ocp_sim_method_num_steps;
ocp_solver = AcadosOcpSolver(ocp);

sim = AcadosSim();
sim.model = model;
sim.model.name = [model_name, '_plant'];

if strcmp(sim_method, 'erk')
	sim.solver_options.integrator_type = 'ERK';
else
	sim.solver_options.integrator_type = 'IRK';
end
sim.solver_options.Tsim = T/ocp_N;
sim.solver_options.num_stages = sim_num_stages;
sim.solver_options.num_steps = sim_num_steps;
sim.solver_options.sens_forw = strcmp(sim_sens_forw, 'true');
sim_solver = AcadosSimSolver(sim);


%% closed loop simulation
N_sim = 200;
x_sim = zeros(nx, N_sim+1);
x_sim(:,1) = x0; % initial state
u_sim = zeros(nu, N_sim);

% set trajectory initialization
%x_traj_init = zeros(nx, ocp_N+1);
%for ii=1:ocp_N x_traj_init(:,ii) = [0; pi; 0; 0]; end
x_traj_init = [linspace(0, 0, ocp_N+1); linspace(pi, 0, ocp_N+1); ...
    linspace(0, 0, ocp_N+1); linspace(0, 0, ocp_N+1)];

u_traj_init = zeros(nu, ocp_N);
pi_traj_init = zeros(nx, ocp_N);



tic;

for ii=1:N_sim

	% set x0
	ocp_solver.set('constr_x0', x_sim(:,ii));

	% set trajectory initialization (if not, set internally using previous solution)
	ocp_solver.set('init_x', x_traj_init);
	ocp_solver.set('init_u', u_traj_init);
	ocp_solver.set('init_pi', pi_traj_init);

	% use ocp_solver.set to modify numerical data for a certain stage
	some_stages = 1:10:ocp_N-1;
	for i = some_stages
        if strcmp(ocp_solver.ocp.cost.cost_type, 'LINEAR_LS')
            ocp_solver.set('cost_Vx', Vx, i);
        end
	end

	% solve OCP
	ocp_solver.solve();

	status = ocp_solver.get('status');
	sqp_iter = ocp_solver.get('sqp_iter');
	time_tot = ocp_solver.get('time_tot');
	time_lin = ocp_solver.get('time_lin');
	time_qp_sol = ocp_solver.get('time_qp_sol');

	fprintf('\nstatus = %d, sqp_iter = %d, time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms])\n',...
		status, sqp_iter, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3);
	if status~=0
		disp('acados ocp solver failed');
		keyboard
	end

	% get solution for initialization of next NLP
	x_traj = ocp_solver.get('x');
	u_traj = ocp_solver.get('u');
	pi_traj = ocp_solver.get('pi');

	% shift trajectory for initialization
	x_traj_init = [x_traj(:,2:end), x_traj(:,end)];
	u_traj_init = [u_traj(:,2:end), u_traj(:,end)];
	pi_traj_init = [pi_traj(:,2:end), pi_traj(:,end)];

	% get solution for sim
	u_sim(:,ii) = ocp_solver.get('u', 0);

	% set initial state of sim
	sim_solver.set('x', x_sim(:,ii));
	% set input in sim
	sim_solver.set('u', u_sim(:,ii));

	% simulate state
	sim_solver.solve();

	% get new state
	x_sim(:,ii+1) = sim_solver.get('xn');

end

avg_time_solve = toc/N_sim


DO_PLOT = 0;
% figures
if DO_PLOT

    for ii=1:N_sim+1
        x_cur = x_sim(:,ii);
    % 	visualize;
    end

    figure;
    subplot(2,1,1);
    plot(0:N_sim, x_sim);
    xlim([0 N_sim]);
    legend('p', 'theta', 'v', 'omega');
    subplot(2,1,2);
    plot(0:N_sim-1, u_sim);
    xlim([0 N_sim]);
    legend('F');


    if is_octave()
        waitforbuttonpress;
    end
end