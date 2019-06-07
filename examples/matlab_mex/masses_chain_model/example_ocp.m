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
%dyn_type = 'explicit';
dyn_type = 'implicit';
%dyn_type = 'discrete';
%sim_method = 'erk';
%sim_method = 'irk';
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
lbu = -10.0*ones(nbu, 1);
ubu =  10.0*ones(nbu, 1);



%% acados ocp model
ocp_model = acados_ocp_model();
% dims
ocp_model.set('T', T);
ocp_model.set('dim_nx', nx);
ocp_model.set('dim_nu', nu);
ocp_model.set('dim_ny', ny);
ocp_model.set('dim_ny_e', ny_e);
if (ng>0)
	ocp_model.set('dim_ng', ng);
	ocp_model.set('dim_ng_e', ng_e);
elseif (nh>0)
	ocp_model.set('dim_nh', nh);
	ocp_model.set('dim_nh_e', nh_e);
else
	ocp_model.set('dim_nbx', nbx);
	ocp_model.set('dim_nbu', nbu);
end
ocp_model.set('dim_np', np);
% symbolics
ocp_model.set('sym_x', model.sym_x);
if isfield(model, 'sym_u')
	ocp_model.set('sym_u', model.sym_u);
end
if isfield(model, 'sym_xdot')
	ocp_model.set('sym_xdot', model.sym_xdot);
end
if (np>0)
	ocp_model.set('sym_p', model.sym_p);
end
% cost
ocp_model.set('cost_type', cost_type);
ocp_model.set('cost_type_e', cost_type);
ocp_model.set('cost_Vu', Vu);
ocp_model.set('cost_Vx', Vx);
ocp_model.set('cost_Vx_e', Vx_e);
ocp_model.set('cost_W', W);
ocp_model.set('cost_W_e', W_e);
ocp_model.set('cost_yr', yr);
ocp_model.set('cost_yr_e', yr_e);
% dynamics
if (strcmp(dyn_type, 'explicit'))
	ocp_model.set('dyn_type', 'explicit');
	ocp_model.set('dyn_expr_f', model.expr_f_expl);
elseif (strcmp(dyn_type, 'implicit'))
	ocp_model.set('dyn_type', 'implicit');
	ocp_model.set('dyn_expr_f', model.expr_f_impl);
else
	ocp_model.set('dyn_type', 'discrete');
	ocp_model.set('dyn_expr_phi', model.expr_phi);
	if (np>0)
		ocp_model.set('dyn_param_phi', 'true');
	end
end
% constraints
ocp_model.set('constr_x0', x0);
if (ng>0)
	ocp_model.set('constr_C', C);
	ocp_model.set('constr_D', D);
	ocp_model.set('constr_lg', lg);
	ocp_model.set('constr_ug', ug);
	ocp_model.set('constr_C_e', C_e);
	ocp_model.set('constr_lg_e', lg_e);
	ocp_model.set('constr_ug_e', ug_e);
elseif (nh>0)
	ocp_model.set('constr_expr_h', model.expr_h);
	ocp_model.set('constr_lh', lh);
	ocp_model.set('constr_uh', uh);
	ocp_model.set('constr_expr_h_e', model.expr_h_e);
	ocp_model.set('constr_lh_e', lh_e);
	ocp_model.set('constr_uh_e', uh_e);
else
	ocp_model.set('constr_Jbx', Jbx);
	ocp_model.set('constr_lbx', lbx);
	ocp_model.set('constr_ubx', ubx);
	ocp_model.set('constr_Jbu', Jbu);
	ocp_model.set('constr_lbu', lbu);
	ocp_model.set('constr_ubu', ubu);
end

%ocp_model.model_struct



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
if (strcmp(dyn_type, 'explicit'))
	ocp_opts.set('sim_method', 'erk');
	ocp_opts.set('sim_method_num_stages', sim_method_num_stages);
	ocp_opts.set('sim_method_num_steps', sim_method_num_steps);
elseif (strcmp(dyn_type, 'implicit'))
	ocp_opts.set('sim_method', 'irk');
	ocp_opts.set('sim_method_num_stages', sim_method_num_stages);
	ocp_opts.set('sim_method_num_steps', sim_method_num_steps);
end

%ocp_opts.opts_struct



%% acados ocp
% create ocp
ocp = acados_ocp(ocp_model, ocp_opts);
%ocp
%ocp.C_ocp
%ocp.C_ocp_ext_fun



% set trajectory initialization
x_traj_init = repmat(model.x_ref, 1, N+1);
u_traj_init = zeros(nu, N);
ocp.set('init_x', x_traj_init);
ocp.set('init_u', u_traj_init);

% set parameter
ocp.set('p', T/N);

% solve
nrep = 1;
tic;
for rep=1:nrep
	ocp.set('init_x', x_traj_init);
	ocp.set('init_u', u_traj_init);
	ocp.solve();
end
time_ext = toc/nrep


% get solution
u = ocp.get('u');
x = ocp.get('x');
%u
%x



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
	fprintf('\niter\tres_g\t\tres_b\t\tres_d\t\tres_m\t\tqp_iter');
	if size(stat,2)>6
		fprintf('\tqp_res_g\tqp_res_b\tqp_res_d\tqp_res_m');
	end
	fprintf('\n');
	for ii=1:size(stat,1)
		fprintf('%d\t%e\t%e\t%e\t%e\t%d', stat(ii,1), stat(ii,2), stat(ii,3), stat(ii,4), stat(ii,5), stat(ii,6));
		if size(stat,2)>6
			fprintf('\t%e\t%e\t%e\t%e', stat(ii,7), stat(ii,8), stat(ii,9), stat(ii,10));
		end
		fprintf('\n');
	end
	fprintf('\n');
else % sqp_rti
	fprintf('\niter\tqp_iter');
	if size(stat,2)>2
		fprintf('\tqp_res_g\tqp_res_b\tqp_res_d\tqp_res_m');
	end
	fprintf('\n');
	for ii=1:size(stat,1)
		fprintf('%d\t%d', stat(ii,1), stat(ii,2));
		if size(stat,2)>2
			fprintf('\t%e\t%e\t%e\t%e', stat(ii,3), stat(ii,4), stat(ii,5), stat(ii,6));
		end
		fprintf('\n');
	end
	fprintf('\n');
end


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


% figures

for ii=1:N
	cur_pos = x(:,ii);
	visualize;
end

if (strcmp(nlp_solver, 'sqp'))
	figure(2);
	plot([0: sqp_iter], log10(stat(:,2)), 'r-x');
	hold on
	plot([0: sqp_iter], log10(stat(:,3)), 'b-x');
	plot([0: sqp_iter], log10(stat(:,4)), 'g-x');
	plot([0: sqp_iter], log10(stat(:,5)), 'k-x');
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
