%% test of native matlab interface
clear all



%% arguments
compile_mex = 'true';
codgen_model = 'true';
param_scheme = 'multiple_shooting_unif_grid';
N = 40;

nlp_solver = 'sqp';
%nlp_solver = 'sqp_rti';
qp_solver = 'partial_condensing_hpipm';
%qp_solver = 'full_condensing_hpipm';
qp_solver_N_pcond = 5;
%sim_method = 'erk';
sim_method = 'irk';
sim_method_num_stages = 4;
sim_method_num_steps = 1;
cost_type = 'linear_ls';



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
nh = 0; % 1
nh_e = 0;
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
%yr = [model.x_ref; zeros(nu, 1)]; % output reference in lagrange term
%yr_e = model.x_ref; % output reference in mayer term

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
Pel_max = 5.0;

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
% TODO power constraint h

% shift
x_end = zeros(nx, 1);
u_end = zeros(nu, 1);

% TODO


%% acados ocp model
ocp_model = acados_ocp_model();
%% dims
ocp_model.set('T', T);
ocp_model.set('nx', nx);
ocp_model.set('nu', nu);
ocp_model.set('ny', ny);
ocp_model.set('ny_e', ny_e);
ocp_model.set('nbx', nbx);
ocp_model.set('nbu', nbu);
%ocp_model.set('nh', nh);
%ocp_model.set('nh_e', nh_e);
ocp_model.set('np', np);
%% symbolics
ocp_model.set('sym_x', model.sym_x);
ocp_model.set('sym_u', model.sym_u);
ocp_model.set('sym_xdot', model.sym_xdot);
ocp_model.set('sym_p', model.sym_p);
%% cost
ocp_model.set('cost_type', cost_type);
ocp_model.set('cost_e_type', cost_type);
ocp_model.set('Vu', Vu);
ocp_model.set('Vx', Vx);
ocp_model.set('Vx_e', Vx_e);
ocp_model.set('W', W);
ocp_model.set('W_e', W_e);
%% dynamics
if (strcmp(sim_method, 'erk'))
	ocp_model.set('dyn_type', 'explicit');
	ocp_model.set('expr_f', model.expr_f_expl);
else % irk
	ocp_model.set('dyn_type', 'implicit');
	ocp_model.set('expr_f', model.expr_f_impl);
end
ocp_model.set('param_f', 'true');
%% constraints
%ocp_model.set('expr_h', model.expr_h);
%ocp_model.set('lh', lh);
%ocp_model.set('uh', uh);
%ocp_model.set('expr_h_e', model.expr_h_e);
%ocp_model.set('lh_e', lh_e);
%ocp_model.set('uh_e', uh_e);
ocp_model.set('Jbx', Jbx);
ocp_model.set('lbx', lbx);
ocp_model.set('ubx', ubx);
ocp_model.set('Jbu', Jbu);
ocp_model.set('lbu', lbu);
ocp_model.set('ubu', ubu);

ocp_model.model_struct



%% acados ocp opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('compile_mex', compile_mex);
ocp_opts.set('codgen_model', codgen_model);
ocp_opts.set('param_scheme', param_scheme);
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('qp_solver', qp_solver);
if (strcmp(qp_solver, 'partial_condensing_hpipm'))
	ocp_opts.set('qp_solver_N_pcond', qp_solver_N_pcond);
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
%ocp.C_ocp_ext_fun



%% solution
% get references
compute_setup;

% set trajectory initialization
x_traj_init = repmat(x0_ref, 1, N+1);
u_traj_init = repmat(u0_ref, 1, N);

ocp.set('x_init', x_traj_init);
ocp.set('u_init', u_traj_init);

% set x0
ocp.set('x0', x0_ref);

% set parameter
nn = 1;
ocp.set('p', wind0_ref(:,nn));

% set reference
ocp.set('yr', y_ref);
ocp.set('yr_e', y_ref);

% solve
ocp.solve();

% get solution
u = ocp.get('u');
x = ocp.get('x');

x(:,1)'
electrical_power = 0.944*97/100*x(1,1)*x(6,1)

status = ocp.get('status');
sqp_iter = ocp.get('sqp_iter');
time_tot = ocp.get('time_tot');
time_lin = ocp.get('time_lin');
time_qp_sol = ocp.get('time_qp_sol');

fprintf('\nstatus = %d, sqp_iter = %d, time_tot = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms])\n', status, sqp_iter, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3);



fprintf('\nsuccess!\n\n');



return;
