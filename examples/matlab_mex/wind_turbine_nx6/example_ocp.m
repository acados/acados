%% test of native matlab interface
clear all



% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	disp('ERROR: env.sh has not been sourced! Before executing this example, run:');
	disp('source env.sh');
	return;
end



%% arguments
compile_mex = 'true';
codgen_model = 'true';
param_scheme = 'multiple_shooting_unif_grid';
N = 40;
nlp_solver = 'sqp';
%nlp_solver = 'sqp_rti';
%nlp_solver_exact_hessian = 'false';
nlp_solver_exact_hessian = 'true';
%regularize_method = 'no_regularize';
%regularize_method = 'project';
regularize_method = 'project_reduc_hess';
%regularize_method = 'mirror';
%regularize_method = 'convexify';
nlp_solver_max_iter = 100;
nlp_solver_ext_qp_res = 1;
qp_solver = 'partial_condensing_hpipm';
%qp_solver = 'full_condensing_hpipm';
qp_solver_cond_N = 5;
qp_solver_cond_ric_alg = 0;
qp_solver_ric_alg = 0;
qp_solver_warm_start = 0;
sim_method = 'erk';
%sim_method = 'irk';
sim_method_num_stages = 4;
if (strcmp(sim_method, 'erk'))
	sim_method_num_steps = 4;
else
	sim_method_num_steps = 1;
end
cost_type = 'linear_ls';
%cost_type = 'nonlinear_ls';



%% create model entries
model = ocp_model_wind_turbine_nx6;



%% dims
Ts = 0.2; % samplig time
T = N*Ts; %8.0; % horizon length time [s]
nx = model.nx; % 8
nu = model.nu; % 2
ny = 4; % number of outputs in lagrange term
ny_e = 2; % number of outputs in mayer term
nbx = 3;
nbu = nu;
ng = 0;
ng_e = 0;
nh = 1;
nh_e = 1;
ns = 1;
ns_e = 1;
nsh = 1;
nsh_e = 1;
np = model.np; % 1

%% cost
% state-to-output matrix in lagrange term
Vx = zeros(ny, nx);
Vx(1, 1) = 1.0;
Vx(2, 5) = 1.0;
% input-to-output matrix in lagrange term
Vu = zeros(ny, nu);
Vu(3, 1) = 1.0;
Vu(4, 2) = 1.0;
% state-to-output matrix in mayer term
Vx_e = zeros(ny_e, nx);
Vx_e(1, 1) = 1.0;
Vx_e(2, 5) = 1.0;
% weight matrix in lagrange term
W = zeros(ny, ny);
W(1, 1) =  1.5114;
W(2, 1) = -0.0649;
W(1, 2) = -0.0649;
W(2, 2) =  0.0180;
W(3, 3) =  0.01;
W(4, 4) =  0.001;
% weight matrix in mayer term
W_e = zeros(ny_e, ny_e); 
W_e(1, 1) =  1.5114;
W_e(2, 1) = -0.0649;
W_e(1, 2) = -0.0649;
W_e(2, 2) =  0.0180;
% output reference in lagrange term
%yr = ... ;
% output reference in mayer term
%yr_e = ... ;
% slacks
Z = 1e2;
Z_e = 1e2;
z = 0e2;
z_e = 0e2;

%% constraints
% constants
dbeta_min = -8.0;
dbeta_max =  8.0;
dM_gen_min = -1.0;
dM_gen_max =  1.0;
OmegaR_min =  6.0/60*2*3.14159265359;
OmegaR_max = 13.0/60*2*3.14159265359;
beta_min =  0.0;
beta_max = 35.0;
M_gen_min = 0.0;
M_gen_max = 5.0;
Pel_min = 0.0;
Pel_max = 5.0; % 5.0

%acados_inf = 1e8;

% state bounds
Jbx = zeros(nbx, nx);
Jbx(1, 1) = 1.0;
Jbx(2, 7) = 1.0;
Jbx(3, 8) = 1.0;
lbx = [OmegaR_min; beta_min; M_gen_min];
ubx = [OmegaR_max; beta_max; M_gen_max];
% input bounds
Jbu = eye(nu);
lbu = [dbeta_min; dM_gen_min];
ubu = [dbeta_max; dM_gen_max];
% nonlinear constraints (power constraint)
lh = Pel_min;
uh = Pel_max;
lh_e = Pel_min;
uh_e = Pel_max;
% soft nonlinear constraints
Jsh = zeros(nh, nsh);
Jsh(1, 1) = 1.0;
Jsh_e = zeros(nh_e, nsh_e);
Jsh_e(1, 1) = 1.0;

% shift
x_end = zeros(nx, 1);
u_end = zeros(nu, 1);



%% acados ocp model
ocp_model = acados_ocp_model();
%% dims
ocp_model.set('T', T);
ocp_model.set('dim_nx', nx);
ocp_model.set('dim_nu', nu);
ocp_model.set('dim_ny', ny);
ocp_model.set('dim_ny_e', ny_e);
ocp_model.set('dim_nbx', nbx);
ocp_model.set('dim_nbu', nbu);
ocp_model.set('dim_nh', nh);
ocp_model.set('dim_nh_e', nh_e);
ocp_model.set('dim_ns', ns);
ocp_model.set('dim_ns_e', ns_e);
ocp_model.set('dim_nsh', nsh);
ocp_model.set('dim_nsh_e', nsh_e);
ocp_model.set('dim_np', np);
%% symbolics
ocp_model.set('sym_x', model.sym_x);
ocp_model.set('sym_u', model.sym_u);
ocp_model.set('sym_xdot', model.sym_xdot);
ocp_model.set('sym_p', model.sym_p);
%% cost
ocp_model.set('cost_type', cost_type);
ocp_model.set('cost_type_e', cost_type);
if (strcmp(cost_type, 'linear_ls'))
	ocp_model.set('cost_Vu', Vu);
	ocp_model.set('cost_Vx', Vx);
	ocp_model.set('cost_Vx_e', Vx_e);
else % nonlinear_ls
	ocp_model.set('cost_expr_y', model.expr_y);
	ocp_model.set('cost_expr_y_e', model.expr_y_e);
end
ocp_model.set('cost_W', W);
ocp_model.set('cost_W_e', W_e);
ocp_model.set('cost_Z', Z);
ocp_model.set('cost_Z_e', Z_e);
ocp_model.set('cost_z', z);
ocp_model.set('cost_z_e', z_e);
%% dynamics
if (strcmp(sim_method, 'erk'))
	ocp_model.set('dyn_type', 'explicit');
	ocp_model.set('dyn_expr_f', model.expr_f_expl);
else % irk
	ocp_model.set('dyn_type', 'implicit');
	ocp_model.set('dyn_expr_f', model.expr_f_impl);
end
% state bounds
ocp_model.set('constr_Jbx', Jbx);
ocp_model.set('constr_lbx', lbx);
ocp_model.set('constr_ubx', ubx);
% input bounds
ocp_model.set('constr_Jbu', Jbu);
ocp_model.set('constr_lbu', lbu);
ocp_model.set('constr_ubu', ubu);
% nonlinear constraints
ocp_model.set('constr_expr_h', model.expr_h);
ocp_model.set('constr_lh', lh);
ocp_model.set('constr_uh', uh);
ocp_model.set('constr_expr_h_e', model.expr_h_e);
ocp_model.set('constr_lh_e', lh_e);
ocp_model.set('constr_uh_e', uh_e);
% soft nonlinear constraints
ocp_model.set('constr_Jsh', Jsh);
ocp_model.set('constr_Jsh_e', Jsh_e);

ocp_model.model_struct



%% acados ocp opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('compile_mex', compile_mex);
ocp_opts.set('codgen_model', codgen_model);
ocp_opts.set('param_scheme', param_scheme);
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('nlp_solver_exact_hessian', nlp_solver_exact_hessian);
ocp_opts.set('regularize_method', regularize_method);
ocp_opts.set('nlp_solver_ext_qp_res', nlp_solver_ext_qp_res);
if (strcmp(nlp_solver, 'sqp'))
	ocp_opts.set('nlp_solver_max_iter', nlp_solver_max_iter);
end
ocp_opts.set('qp_solver', qp_solver);
if (strcmp(qp_solver, 'partial_condensing_hpipm'))
	ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
	ocp_opts.set('qp_solver_cond_ric_alg', qp_solver_cond_ric_alg);
	ocp_opts.set('qp_solver_ric_alg', qp_solver_ric_alg);
	ocp_opts.set('qp_solver_warm_start', qp_solver_warm_start);
end
ocp_opts.set('sim_method', sim_method);
ocp_opts.set('sim_method_num_stages', sim_method_num_stages);
ocp_opts.set('sim_method_num_steps', sim_method_num_steps);

ocp_opts.opts_struct



%% acados ocp
% create ocp
ocp = acados_ocp(ocp_model, ocp_opts);
%ocp
%ocp.C_ocp
ocp.C_ocp_ext_fun



%% solution
% get references
compute_setup;

% set trajectory initialization
x_traj_init = repmat(x0_ref, 1, N+1);
u_traj_init = repmat(u0_ref, 1, N);

tic

ocp.set('init_x', x_traj_init);
ocp.set('init_u', u_traj_init);

% set x0
ocp.set('constr_x0', x0_ref);

% set parameter
nn = 1;
ocp.set('p', wind0_ref(:,nn));

% set reference
ocp.set('cost_yr', y_ref(:,nn));
ocp.set('cost_yr_e', y_ref(:,nn));

% solve
disp('before solve')
ocp.solve();
disp('after solve')

% get solution
u = ocp.get('u');
x = ocp.get('x');

time_ext = toc;

x(:,1)'
u(:,1)'
%electrical_power = 0.944*97/100*x(1,1)*x(6,1)
electrical_power = 0.944*97/100*x(1,:).*x(6,:)


% statistics

status = ocp.get('status');
sqp_iter = ocp.get('sqp_iter');
time_tot = ocp.get('time_tot');
time_lin = ocp.get('time_lin');
time_reg = ocp.get('time_reg');
time_qp_sol = ocp.get('time_qp_sol');

fprintf('\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n', status, sqp_iter, time_ext*1e3, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3);

stat = ocp.get('stat');
if (strcmp(nlp_solver, 'sqp'))
	fprintf('\niter\tres_g\t\tres_b\t\tres_d\t\tres_m\t\tqp_stat\tqp_iter');
	if size(stat,2)>7
		fprintf('\tqp_res_g\tqp_res_b\tqp_res_d\tqp_res_m');
	end
	fprintf('\n');
	for ii=1:size(stat,1)
		fprintf('%d\t%e\t%e\t%e\t%e\t%d\t%d', stat(ii,1), stat(ii,2), stat(ii,3), stat(ii,4), stat(ii,5), stat(ii,6), stat(ii,7));
		if size(stat,2)>7
			fprintf('\t%e\t%e\t%e\t%e', stat(ii,8), stat(ii,9), stat(ii,10), stat(ii,11));
		end
		fprintf('\n');
	end
	fprintf('\n');
else % sqp_rti
	fprintf('\niter\tqp_stat\tqp_iter');
	if size(stat,2)>3
		fprintf('\tqp_res_g\tqp_res_b\tqp_res_d\tqp_res_m');
	end
	fprintf('\n');
	for ii=1:size(stat,1)
		fprintf('%d\t%d\t%d', stat(ii,1), stat(ii,2), stat(ii,3));
		if size(stat,2)>3
			fprintf('\t%e\t%e\t%e\t%e', stat(ii,4), stat(ii,5), stat(ii,6), stat(ii,7));
		end
		fprintf('\n');
	end
	fprintf('\n');
end


% figures

figure(1);
subplot(3,1,1);
plot(0:N, x);
xlim([0 N]);
ylabel('states');
%legend('p', 'theta', 'v', 'omega');
subplot(3,1,2);
plot(0:N-1, u);
xlim([0 N]);
ylabel('inputs');
%legend('F');
subplot(3,1,3);
plot(0:N, electrical_power);
hold on
plot([0 N], [Pel_max Pel_max]);
hold off
xlim([0 N]);
ylim([4.5 6.0]);
ylabel('electrical power');
%legend('F');

if (strcmp(nlp_solver, 'sqp'))
	figure(2);
	plot([0: size(stat,1)-1], log10(stat(:,2)), 'r-x');
	hold on
	plot([0: size(stat,1)-1], log10(stat(:,3)), 'b-x');
	plot([0: size(stat,1)-1], log10(stat(:,4)), 'g-x');
	plot([0: size(stat,1)-1], log10(stat(:,5)), 'k-x');
	hold off
	xlabel('iter')
	ylabel('res')
end



if status==0
	fprintf('\nsuccess!\n\n');
else
	fprintf('\nsolution failed!\n\n');
end


waitforbuttonpress;


return;
