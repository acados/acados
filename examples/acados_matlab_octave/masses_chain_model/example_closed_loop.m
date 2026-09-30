%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.




%% example of closed loop simulation
clear all

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end

%% handy arguments
compile_interface = 'auto';
% simulation
sim_method = 'irk';
sim_sens_forw = 'false';
sim_num_stages = 4;
sim_num_steps = 4;
% ocp
ocp_N = 40;
%ocp_nlp_solver = 'sqp';
ocp_nlp_solver = 'sqp_rti';
ocp_nlp_solver_exact_hessian = 'false';
%ocp_nlp_solver_exact_hessian = 'true';
regularize_method = 'no_regularize';
%regularize_method = 'project';
%regularize_method = 'project_reduc_hess';
%regularize_method = 'mirror';
%regularize_method = 'convexify';
ocp_nlp_solver_max_iter = 100;
ocp_nlp_solver_ext_qp_res = 1;
ocp_nlp_solver_warm_start_first_qp = 1;
ocp_qp_solver = 'partial_condensing_hpipm';
%ocp_qp_solver = 'full_condensing_hpipm';
%ocp_qp_solver = 'full_condensing_qpoases';
%ocp_qp_solver = 'partial_condensing_osqp';
ocp_qp_solver_cond_N = 5;
%ocp_qp_solver_cond_N = ocp_N;
ocp_qp_solver_cond_ric_alg = 0;
ocp_qp_solver_ric_alg = 0;
ocp_qp_solver_warm_start = 1;
ocp_qp_solver_max_iter = 50;
%ocp_sim_method = 'erk';
ocp_sim_method = 'irk';
ocp_sim_method_num_stages = 4;
ocp_sim_method_num_steps = 2;
ocp_cost_type = 'linear_ls';


%% create model entries
nfm = 4;    % number of free masses
nm = nfm+1; % number of masses
model = masses_chain_model(nfm);
model_name = 'masses_chain_closed_loop';
wall = -0.01;


% dims
T = 8.0; % horizon length time
nx = model.nx; % 6*nfm
nu = model.nu; % 3
ny = nu+nx; % number of outputs in lagrange term
ny_e = nx; % number of outputs in mayer term
nbx = nfm;
nbu = nu;
ng = 0;
nh = 0;
nh_e = 0;

% cost
Vx = zeros(ny, nx); for ii=1:nx Vx(ii,ii)=1.0; end % state-to-output matrix in lagrange term
Vu = zeros(ny, nu); for ii=1:nu Vu(nx+ii,ii)=1.0; end % input-to-output matrix in lagrange term
Vx_e = zeros(ny_e, nx); for ii=1:nx Vx_e(ii,ii)=1.0; end % state-to-output matrix in mayer term
W = 10.0*eye(ny); for ii=1:nu W(nx+ii,nx+ii)=1e-2; end % weight matrix in lagrange term
W_e = 10.0*eye(ny_e); % weight matrix in mayer term
yr = [model.x_ref; zeros(nu, 1)]; % output reference in lagrange term
yr_e = model.x_ref; % output reference in mayer term

% constraints
x0 = model.x0;
%x0 = model.x_ref;
Jbx = zeros(nbx, nx); for ii=1:nbx Jbx(ii,2+6*(ii-1))=1.0; end
lbx = wall*ones(nbx, 1);
ubx = 1e+4*ones(nbx, 1);
Jbu = zeros(nbu, nu); for ii=1:nbu Jbu(ii,ii)=1.0; end
lbu = -1.0*ones(nbu, 1);
ubu =  1.0*ones(nbu, 1);


%% OCP and simulation formulations
ocp = AcadosOcp();
ocp.model.name = model_name;
ocp.model.x = model.sym_x;
ocp.model.u = model.sym_u;
ocp.model.xdot = model.sym_xdot;
ocp.model.f_impl_expr = model.expr_f_impl;
ocp.model.f_expl_expr = model.expr_f_expl;

ocp.cost.cost_type_0 = 'LINEAR_LS';
ocp.cost.cost_type = 'LINEAR_LS';
ocp.cost.cost_type_e = 'LINEAR_LS';
ocp.cost.Vu_0 = Vu;
ocp.cost.Vx_0 = Vx;
ocp.cost.W_0 = W;
ocp.cost.yref_0 = yr;
ocp.cost.Vu = Vu;
ocp.cost.Vx = Vx;
ocp.cost.W = W;
ocp.cost.yref = yr;
ocp.cost.Vx_e = Vx_e;
ocp.cost.W_e = W_e;
ocp.cost.yref_e = yr_e;

ocp.constraints.x0 = x0;
ocp.constraints.idxbx = (1:6:nx)';
ocp.constraints.lbx = lbx;
ocp.constraints.ubx = ubx;
ocp.constraints.idxbu = (0:nu-1)';
ocp.constraints.lbu = lbu;
ocp.constraints.ubu = ubu;

ocp.solver_options.N_horizon = ocp_N;
ocp.solver_options.tf = T;
ocp.solver_options.nlp_solver_type = upper(ocp_nlp_solver);
if strcmp(ocp_nlp_solver_exact_hessian, 'true')
	ocp.solver_options.hessian_approx = 'EXACT';
else
	ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
end
ocp.solver_options.regularize_method = upper(regularize_method);
ocp.solver_options.nlp_solver_ext_qp_res = ocp_nlp_solver_ext_qp_res;
ocp.solver_options.nlp_solver_warm_start_first_qp = ocp_nlp_solver_warm_start_first_qp;
ocp.solver_options.qp_solver = upper(ocp_qp_solver);
ocp.solver_options.integrator_type = upper(ocp_sim_method);
ocp.solver_options.qp_solver_iter_max = ocp_qp_solver_max_iter;
ocp.solver_options.qp_solver_warm_start = ocp_qp_solver_warm_start;
if contains(ocp_qp_solver, 'partial_condensing')
	ocp.solver_options.qp_solver_cond_N = ocp_qp_solver_cond_N;
end
if strcmp(ocp_qp_solver, 'partial_condensing_hpipm')
	ocp.solver_options.qp_solver_cond_ric_alg = ocp_qp_solver_cond_ric_alg;
	ocp.solver_options.qp_solver_ric_alg = ocp_qp_solver_ric_alg;
end
ocp.solver_options.sim_method_num_stages = ocp_sim_method_num_stages;
ocp.solver_options.sim_method_num_steps = ocp_sim_method_num_steps;

ocp_solver = AcadosOcpSolver(ocp);

sim = AcadosSim();
sim.model.name = [model_name, '_plant'];
sim.model.x = model.sym_x;
sim.model.u = model.sym_u;
sim.model.xdot = model.sym_xdot;
sim.model.f_impl_expr = model.expr_f_impl;
sim.model.f_expl_expr = model.expr_f_expl;
sim.solver_options.Tsim = T/ocp_N;
sim.solver_options.integrator_type = upper(sim_method);
sim.solver_options.num_stages = sim_num_stages;
sim.solver_options.num_steps = sim_num_steps;
sim.solver_options.sens_forw = strcmp(sim_sens_forw, 'true');
sim_solver = AcadosSimSolver(sim);


%% closed loop simulation
n_sim = 50;
x_sim = zeros(nx, n_sim+1);
x_sim(:,1) = x0; % initial state
u_sim = zeros(nu, n_sim);

x_traj_init = repmat(model.x_ref, 1, ocp_N+1);
u_traj_init = zeros(nu, ocp_N);
pi_traj_init = zeros(nx, ocp_N);

%ocp_solver.set('init_x', x_traj_init);
%ocp_solver.set('init_u', u_traj_init);
%ocp_solver.set('init_pi', pi_traj_init);

tic;

for ii=1:n_sim

	% set x0
	ocp_solver.set('constr_x0', x_sim(:,ii));

	% set trajectory initialization (if not, set internally using previous solution)
	ocp_solver.set('init_x', x_traj_init);
	ocp_solver.set('init_u', u_traj_init);
	ocp_solver.set('init_pi', pi_traj_init);

	% solve OCP
	ocp_solver.set('rti_phase', 1);
	ocp_solver.solve();
	ocp_solver.set('rti_phase', 2);
	ocp_solver.solve();

	if 1
		status = ocp_solver.get('status');
		sqp_iter = ocp_solver.get('sqp_iter');
		time_tot = ocp_solver.get('time_tot');
		time_lin = ocp_solver.get('time_lin');
		time_reg = ocp_solver.get('time_reg');
		time_qp_sol = ocp_solver.get('time_qp_sol');
		time_qp_solver_call = ocp_solver.get('time_qp_solver_call');
		qp_iter = ocp_solver.get('qp_iter_all');

		fprintf('\nstatus = %d, sqp_iter = %d, time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms] (time_qp_solver_call = %f [ms]), time_reg = %f [ms])\n', status, sqp_iter, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_qp_solver_call*1e3, time_reg*1e3);
%		fprintf('%e %d\n', time_qp_solver_call, qp_iter);

%		ocp_solver.print('stat');
	end

	% get solution
	x_traj = ocp_solver.get('x');
	u_traj = ocp_solver.get('u');
	pi_traj = ocp_solver.get('pi');

	% shift trajectory for initialization
	x_traj_init = [x_traj(:,2:end), x_traj(:,end)];
	u_traj_init = [u_traj(:,2:end), u_traj(:,end)];
	pi_traj_init = [pi_traj(:,2:end), pi_traj(:,end)];

	% get solution for sim
	u_sim(:,ii) = ocp_solver.get('u', 0);

	% overwrite control to perturb the system
%	if(ii<=5)
%		u_sim(:,ii) = [-1; 1; 1];
%	end

	% set initial state of sim
	sim_solver.set('x', x_sim(:,ii));
	% set input in sim
	sim_solver.set('u', u_sim(:,ii));

	% simulate state
	sim_solver.solve();

	% get new state
	x_sim(:,ii+1) = sim_solver.get('xn');

end

avg_time_solve = toc/n_sim


u_sim;
x_sim;

% print solution
for ii=1:n_sim+1
	cur_pos = x_sim(:,ii);
	visualize;
end



status = ocp_solver.get('status');

if status==0
	fprintf('\nsuccess!\n\n');
else
	fprintf('\nsolution failed!\n\n');
end


if is_octave()
    waitforbuttonpress;
end
