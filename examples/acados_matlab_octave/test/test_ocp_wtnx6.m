%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

%% test of native matlab interface

addpath('../wind_turbine_nx6/');

% qpOASES too slow to test on CI
for itest = 1:3

%% arguments
% simulation
sim_method = 'irk';
sim_sens_forw = 'false';
sim_num_stages = 4;
sim_num_steps = 1;
% ocp
ocp_N = 40;
ocp_nlp_solver = 'sqp';
%ocp_nlp_solver = 'sqp_rti';
ocp_nlp_solver_exact_hessian = 'false';
%ocp_nlp_solver_exact_hessian = 'true';
regularize_method = 'no_regularize';
ocp_nlp_solver_max_iter = 50;
ocp_nlp_solver_tol_stat = 1e-8;
ocp_nlp_solver_tol_eq   = 1e-8;
ocp_nlp_solver_tol_ineq = 1e-8;
ocp_nlp_solver_tol_comp = 1e-8;
ocp_nlp_solver_ext_qp_res = 1;
switch itest
case 1
    ocp_qp_solver = 'partial_condensing_hpipm';
case 2
    ocp_qp_solver = 'full_condensing_hpipm';
case 3
    ocp_qp_solver = 'full_condensing_daqp';
case 4
    ocp_qp_solver = 'full_condensing_qpoases';
end
fprintf(['\n\nrunning with qp solver ', ocp_qp_solver, '\n'])
ocp_qp_solver_cond_N = 5;
ocp_qp_solver_cond_ric_alg = 0;
ocp_qp_solver_ric_alg = 0;
ocp_qp_solver_warm_start = 2;
%ocp_sim_method = 'erk';
ocp_sim_method = 'irk';
ocp_sim_method_num_stages = 4 * ones(ocp_N, 1); % scalar or vector of size ocp_N;
ocp_sim_method_num_steps = 1 * ones(ocp_N, 1); % scalar or vector of size ocp_N;
ocp_sim_method_newton_iter = 3; % * ones(ocp_N, 1); % scalar or vector of size ocp_N;
%cost_type = 'linear_ls';
cost_type = 'nonlinear_ls';

% get references
compute_setup;


%% create model entries
model = ocp_model_wind_turbine_nx6;



%% dims
Ts = 0.2; % samplig time
T = ocp_N*Ts; %8.0; % horizon length time [s]
nx = length(model.x); % 8
nu = length(model.u); % 2
ny = 4; % number of outputs in lagrange term
ny_e = 2; % number of outputs in mayer term
ns = 2;
ns_e = 2;
np = length(model.p); % 1

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
z = 0e2*ones(ns,1);
z_e = 0e2*ones(ns_e,1);

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

lbx = [OmegaR_min; beta_min; M_gen_min];
ubx = [OmegaR_max; beta_max; M_gen_max];
% input bounds
lbu = [dbeta_min; dM_gen_min];
ubu = [dbeta_max; dM_gen_max];
% nonlinear constraints (power constraint)
lh = Pel_min;
uh = Pel_max;
lh_e = Pel_min;
uh_e = Pel_max;
%% OCP formulation
ocp = AcadosOcp();
ocp.model = model;

model.cost_y_expr_0 = [model.x(1); model.x(5); model.u];
model.cost_y_expr = [model.x(1); model.x(5); model.u];
model.cost_y_expr_e = [model.x(1); model.x(5)];
ocp.cost.cost_type_0 = upper(cost_type);
ocp.cost.cost_type = upper(cost_type);
ocp.cost.cost_type_e = upper(cost_type);
ocp.cost.W_0 = W;
ocp.cost.W = W;
ocp.cost.W_e = W_e;
ocp.cost.yref_0 = zeros(ny, 1);
ocp.cost.yref = zeros(ny, 1);
ocp.cost.yref_e = zeros(ny_e, 1);
if strcmp(cost_type, 'linear_ls')
    ocp.cost.Vx_0 = Vx;
    ocp.cost.Vu_0 = Vu;
    ocp.cost.Vx = Vx;
    ocp.cost.Vu = Vu;
    ocp.cost.Vx_e = Vx_e;
end
ocp.cost.Zl = 1e2 * ones(ns, 1);
ocp.cost.Zu = 1e2 * ones(ns, 1);
ocp.cost.zl = z;
ocp.cost.zu = z;
ocp.cost.Zl_e = 0 * ones(ns_e, 1);;
ocp.cost.Zu_e = 0 * ones(ns_e, 1);;
ocp.cost.zl_e = z_e;
ocp.cost.zu_e = z_e;

ocp.constraints.idxbx = [0; 6; 7];
ocp.constraints.lbx = lbx;
ocp.constraints.ubx = ubx;
ocp.constraints.idxbx_e = [0; 6; 7];
ocp.constraints.lbx_e = lbx;
ocp.constraints.ubx_e = ubx;
ocp.constraints.idxbu = (0:nu-1)';
ocp.constraints.lbu = lbu;
ocp.constraints.ubu = ubu;
ocp.constraints.lh = lh;
ocp.constraints.uh = uh;
ocp.constraints.lh_e = lh_e;
ocp.constraints.uh_e = uh_e;
ocp.constraints.idxsbx = 0;
ocp.constraints.idxsbx_e = 0;
ocp.constraints.idxsh = 0;
ocp.constraints.idxsh_e = 0;
ocp.constraints.x0 = x0_ref;
ocp.parameter_values = wind0_ref(:,1);

ocp.solver_options.N_horizon = ocp_N;
ocp.solver_options.tf = T;
ocp.solver_options.nlp_solver_type = upper(ocp_nlp_solver);
if strcmp(ocp_nlp_solver_exact_hessian, 'false')
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
else
    ocp.solver_options.hessian_approx = 'EXACT';
end
ocp.solver_options.nlp_solver_ext_qp_res = ocp_nlp_solver_ext_qp_res;
ocp.solver_options.nlp_solver_max_iter = ocp_nlp_solver_max_iter;
ocp.solver_options.nlp_solver_tol_stat = ocp_nlp_solver_tol_stat;
ocp.solver_options.nlp_solver_tol_eq = ocp_nlp_solver_tol_eq;
ocp.solver_options.nlp_solver_tol_ineq = ocp_nlp_solver_tol_ineq;
ocp.solver_options.nlp_solver_tol_comp = ocp_nlp_solver_tol_comp;
ocp.solver_options.qp_solver = upper(ocp_qp_solver);
ocp.solver_options.qp_solver_iter_max = 500;
if strcmp(ocp_qp_solver, 'partial_condensing_hpipm')
    ocp.solver_options.qp_solver_cond_N = ocp_qp_solver_cond_N;
    ocp.solver_options.qp_solver_cond_ric_alg = ocp_qp_solver_cond_ric_alg;
    ocp.solver_options.qp_solver_ric_alg = ocp_qp_solver_ric_alg;
    ocp.solver_options.qp_solver_warm_start = ocp_qp_solver_warm_start;
end
ocp.solver_options.integrator_type = upper(ocp_sim_method);
ocp.solver_options.sim_method_num_stages = ocp_sim_method_num_stages;
ocp.solver_options.sim_method_num_steps = ocp_sim_method_num_steps;
ocp.solver_options.sim_method_newton_iter = ocp_sim_method_newton_iter;
ocp.solver_options.regularize_method = 'NO_REGULARIZE';

ocp_solver = AcadosOcpSolver(ocp);
%ocp
%ocp_solver.C_ocp

%% Plant integrator
sim = AcadosSim();
sim.model = model;
sim.solver_options.Tsim = T/ocp_N;
sim.solver_options.integrator_type = upper(sim_method);
sim.solver_options.num_stages = sim_num_stages;
sim.solver_options.num_steps = sim_num_steps;
sim.solver_options.sens_forw = strcmp(sim_sens_forw, 'true');
sim.parameter_values = zeros(np, 1);
sim_solver = AcadosSimSolver(sim);


%% closed loop simulation
n_sim = 100;
n_sim_max = length(wind0_ref) - ocp_N;
if n_sim>n_sim_max
    n_sim = n_sim_max;
end
x_sim = zeros(nx, n_sim+1);
x_sim(:,1) = x0_ref; % initial state
u_sim = zeros(nu, n_sim);

sqp_iter_sim = zeros(n_sim,1);
time_ext = zeros(n_sim, 1);
time_tot = zeros(n_sim, 1);
time_lin = zeros(n_sim, 1);
time_qp_sol = zeros(n_sim, 1);

% set trajectory initialization
x_traj_init = repmat(x0_ref, 1, ocp_N+1);
u_traj_init = repmat(u0_ref, 1, ocp_N);
pi_traj_init = zeros(nx, ocp_N);

for ii=1:n_sim

%    fprintf('\nsimulation step %d\n', ii);

    tic

    % set x0
    ocp_solver.set('constr_x0', x_sim(:,ii));
    % set parameter
    for jj=0:ocp_N-1
        ocp_solver.set('p', wind0_ref(:,ii+jj), jj);
    end

    % set reference (different at each stage)
    for jj=0:ocp_N-1
        ocp_solver.set('cost_y_ref', y_ref(:,ii+jj), jj);
    end
    ocp_solver.set('cost_y_ref', y_ref(1:ny_e,ii+ocp_N), ocp_N);

    % set trajectory initialization (if not, set internally using previous solution)
    ocp_solver.set('init_x', x_traj_init);
    ocp_solver.set('init_u', u_traj_init);
    ocp_solver.set('init_pi', pi_traj_init);

    % solve
    ocp_solver.solve();

    % get solution
    x = ocp_solver.get('x');
    u = ocp_solver.get('u');
    pi = ocp_solver.get('pi');

    % store first input
    u_sim(:,ii) = ocp_solver.get('u', 0);

    % set initial state of sim
    sim_solver.set('x', x_sim(:,ii));
    % set input in sim
    sim_solver.set('u', u_sim(:,ii));
    % set parameter
    sim_solver.set('p', wind0_ref(:,ii));

    % simulate state
    sim_solver.solve();

    % get new state
    x_sim(:,ii+1) = sim_solver.get('xn');

    % shift trajectory for initialization
%    x_traj_init = [x(:,2:ocp_N+1), zeros(nx, 1)];
    x_traj_init = [x(:,2:ocp_N+1), x(:,ocp_N+1)];
%    x_traj_init = [x(:,2:ocp_N+1), sim_solver.get('xn')];
%    u_traj_init = [u(:,2:ocp_N), zeros(nu, 1)];
    u_traj_init = [u(:,2:ocp_N), u(:,ocp_N)];
    pi_traj_init = [pi(:,2:ocp_N), pi(:,ocp_N)];

    time_ext(ii) = toc;

    electrical_power = 0.944*97/100*x(1,1)*x(6,1);

    status = ocp_solver.get('status');
    sqp_iter = ocp_solver.get('sqp_iter');
    time_tot(ii) = ocp_solver.get('time_tot');
    time_lin(ii) = ocp_solver.get('time_lin');
    time_qp_sol(ii) = ocp_solver.get('time_qp_sol');
    sqp_stats = ocp_solver.get('stat');
    qp_iter = sqp_stats(:,7);

    sqp_iter_sim(ii) = sqp_iter;
    if status ~= 0
        ocp_solver.print()
        error(['ocp_nlp solver returned status ', num2str(status), '!= 0 in simulation instance ', num2str(ii)]);
    end

    fprintf('\nstatus = %d, sqp_iter = %d, qp_iter = %d, time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms]), Pel = %f',...
            status, sqp_iter, sum(qp_iter), time_ext(ii)*1e3, time_tot(ii)*1e3, time_lin(ii)*1e3, time_qp_sol(ii)*1e3, electrical_power);

    if 0
        ocp_solver.print('stat')
    end

end

% get slack values
for i = 0:ocp_N-1
    sl = ocp_solver.get('sl', i);
    su = ocp_solver.get('su', i);
    % test setters
    ocp_solver.set('sl', sl, i);
    ocp_solver.set('su', su, i);
end
sl = ocp_solver.get('sl', ocp_N);
su = ocp_solver.get('su', ocp_N);



electrical_power = 0.944*97/100*x_sim(1,:).*x_sim(6,:);

x_sim_ref = [   1.263425730522397
   0.007562725557589
  76.028289356099236
   0.007188510774546
   6.949049234224142
   3.892712459979240
   6.302629591941585
   3.882220255648666];

err_vs_ref = x_sim_ref - x_sim(:,end);

fprintf('\nmedian computation times: time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms])', ...
        median(time_ext)*1e3, median(time_tot)*1e3, median(time_lin)*1e3, median(time_qp_sol)*1e3)

if status~=0
    error('test_ocp_wtnx6: solution failed!');
elseif max(abs(err_vs_ref)) > 1e-14
    error('test_ocp_wtnx6: to high deviation from known result!');
elseif sqp_iter > 2
    error('test_ocp_wtnx6: sqp_iter > 2, this problem is typically solved within less iterations!');
else
    fprintf('\ntest_ocp_wtnx6: success!\n');
end

% figures
if 0
    figure;
    subplot(3,1,1);
    plot(0:n_sim, x_sim);
    xlim([0 n_sim]);
    ylabel('states');
    %legend('p', 'theta', 'v', 'omega');
    subplot(3,1,2);
    plot(0:n_sim-1, u_sim);
    xlim([0 n_sim]);
    ylabel('inputs');
    %legend('F');
    subplot(3,1,3);
    plot(0:n_sim, electrical_power);
    hold on
    plot([0 n_sim], [Pel_max Pel_max]);
    hold off
    xlim([0 n_sim]);
    ylim([4.0 6.0]);
    ylabel('electrical power');
    %legend('F');

    figure;
    plot(1:n_sim, sqp_iter_sim, 'rx');
    hold on
    plot([1 n_sim], [ocp_nlp_solver_max_iter ocp_nlp_solver_max_iter]);
    hold off
    ylim([0 ocp_nlp_solver_max_iter+1])
    ylabel('sqp iterations')
    xlabel('sqp calls')
    if is_octave()
        waitforbuttonpress;
    end
end

end

% remove temporary created files
delete('y_ref.mat')
delete('y_e_ref.mat')
delete('wind0_ref.mat')
delete('windN_ref.mat')
