%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.




clear all


% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end


%% arguments
compile_interface = 'auto';
gnsf_detect_struct = 'true';
model_name = 'masses_chain';

N = 40;
nlp_solver = 'sqp';
%nlp_solver = 'sqp_rti';
nlp_solver_exact_hessian = 'false';
%nlp_solver_exact_hessian = 'true';
regularize_method = 'no_regularize';
%regularize_method = 'project';
%regularize_method = 'project_reduc_hess';
%regularize_method = 'mirror';
%regularize_method = 'convexify';
nlp_solver_max_iter = 100;
nlp_solver_ext_qp_res = 1;
nlp_solver_warm_start_first_qp = 0;
qp_solver = 'partial_condensing_hpipm';
%qp_solver = 'full_condensing_hpipm';
%qp_solver = 'full_condensing_qpoases';
%qp_solver = 'partial_condensing_osqp';
qp_solver_cond_N = 5;
qp_solver_cond_ric_alg = 0;
qp_solver_ric_alg = 0;
qp_solver_warm_start = 0;
qp_solver_max_iter = 100;
%dyn_type = 'explicit';
dyn_type = 'implicit';
%dyn_type = 'discrete';
%sim_method = 'erk';
sim_method = 'irk';
%sim_method = 'irk_gnsf';
sim_method_num_stages = 4;
sim_method_num_steps = 2;
cost_type = 'linear_ls';



%% create model entries
nfm = 4;    % number of free masses
nm = nfm+1; % number of masses
model = masses_chain_model(nfm);
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
ng_e = 0;
nh = 0;
nh_e = 0;
np = model.np;

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
Jbx = zeros(nbx, nx); for ii=1:nbx Jbx(ii,2+6*(ii-1))=1.0; end
lbx = wall*ones(nbx, 1);
ubx = 1e+4*ones(nbx, 1);
Jbu = zeros(nbu, nu); for ii=1:nbu Jbu(ii,ii)=1.0; end
lbu = -1.0*ones(nbu, 1);
ubu =  1.0*ones(nbu, 1);


%% acados OCP
ocp = AcadosOcp();
ocp.model.name = model_name;
ocp.model.x = model.sym_x;
ocp.model.u = model.sym_u;
ocp.model.xdot = model.sym_xdot;
ocp.model.p = model.sym_p;

switch dyn_type
	case 'explicit'
		ocp.model.f_expl_expr = model.expr_f_expl;
		ocp.solver_options.integrator_type = 'ERK';
	case 'implicit'
		ocp.model.f_impl_expr = model.expr_f_impl;
		ocp.solver_options.integrator_type = 'IRK';
	otherwise
		ocp.model.disc_dyn_expr = model.expr_phi;
		ocp.solver_options.integrator_type = 'DISCRETE';
end

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

ocp.solver_options.N_horizon = N;
ocp.solver_options.tf = T;
ocp.solver_options.nlp_solver_type = upper(nlp_solver);
if strcmp(nlp_solver_exact_hessian, 'true')
	ocp.solver_options.hessian_approx = 'EXACT';
else
	ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
end
ocp.solver_options.regularize_method = upper(regularize_method);
ocp.solver_options.nlp_solver_ext_qp_res = nlp_solver_ext_qp_res;
ocp.solver_options.nlp_solver_warm_start_first_qp = nlp_solver_warm_start_first_qp;
ocp.solver_options.nlp_solver_max_iter = nlp_solver_max_iter;
ocp.solver_options.qp_solver = upper(qp_solver);
ocp.solver_options.qp_solver_iter_max = qp_solver_max_iter;
ocp.solver_options.qp_solver_warm_start = qp_solver_warm_start;
if contains(qp_solver, 'partial_condensing')
	ocp.solver_options.qp_solver_cond_N = qp_solver_cond_N;
end
if strcmp(qp_solver, 'partial_condensing_hpipm')
	ocp.solver_options.qp_solver_cond_ric_alg = qp_solver_cond_ric_alg;
	ocp.solver_options.qp_solver_ric_alg = qp_solver_ric_alg;
end
ocp.solver_options.sim_method_num_stages = sim_method_num_stages;
ocp.solver_options.sim_method_num_steps = sim_method_num_steps;
ocp.solver_options.compile_interface = [];

if strcmp(sim_method, 'irk_gnsf')
	ocp.solver_options.integrator_type = 'GNSF';
end

ocp.parameter_values = T/N;
ocp_solver = AcadosOcpSolver(ocp);

% set trajectory initialization
x_traj_init = repmat(model.x_ref, 1, N+1);
u_traj_init = zeros(nu, N);
ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);

% set parameter
ocp_solver.set('p', T/N, 0, N+1);

% solve
nrep = 1;
tic;
for rep=1:nrep
	ocp_solver.set('init_x', x_traj_init);
	ocp_solver.set('init_u', u_traj_init);
	ocp_solver.solve();
end
time_ext = toc/nrep


% get solution
u = ocp_solver.get('u');
x = ocp_solver.get('x');


% statistics
status = ocp_solver.get('status');
sqp_iter = ocp_solver.get('sqp_iter');
time_tot = ocp_solver.get('time_tot');
time_lin = ocp_solver.get('time_lin');
time_reg = ocp_solver.get('time_reg');
time_qp_sol = ocp_solver.get('time_qp_sol');
time_qp_solver_call = ocp_solver.get('time_qp_solver_call');

fprintf('\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms] (time_qp_solver_call = %f [ms]), time_reg = %f [ms])\n', status, sqp_iter, time_ext*1e3, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_qp_solver_call*1e3, time_reg*1e3);

ocp_solver.print('stat');


%% figures
% plot result
%figure()
%subplot(2, 1, 1)
%plot(0:N, x);
%title('closed loop simulation')
%ylabel('x')
%subplot(2, 1, 2)
%plot(1:N, u);
%ylabel('u')
%xlabel('sample')

for ii=1:N
	cur_pos = x(:,ii);
	visualize;
end

stat = ocp_solver.get('stat');
if (strcmp(nlp_solver, 'sqp'))
	figure();
	plot(0: size(stat,1)-1, log10(stat(:,2)), 'r-x');
	hold on
	plot(0: size(stat,1)-1, log10(stat(:,3)), 'b-x');
	plot(0: size(stat,1)-1, log10(stat(:,4)), 'g-x');
	plot(0: size(stat,1)-1, log10(stat(:,5)), 'k-x');
	hold off
	xlabel('iter')
	ylabel('residuals')
    legend('res stat', 'res eq', 'res ineq', 'res compl');
end


if status==0
	fprintf('\nsuccess!\n\n');
else
	fprintf('\nsolution failed!\n\n');
end


if is_octave()
    waitforbuttonpress;
end
