%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

clear all; clc;

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
    error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end

%% arguments
compile_interface = 'auto'; %'auto';
gnsf_detect_struct = 'true';

% discretization
N = 20;
T = 1; % horizon length time
h = T/N;

nlp_solver = 'sqp';
%nlp_solver = 'sqp_rti';
%nlp_solver_exact_hessian = 'false';
nlp_solver_exact_hessian = 'true';
%regularize_method = 'no_regularize';
%regularize_method = 'project';
regularize_method = 'project_reduc_hess';
%regularize_method = 'mirror';
%regularize_method = 'convexify';
nlp_solver_max_iter = 100; %100;
nlp_solver_tol_stat = 1e-8;
nlp_solver_tol_eq   = 1e-8;
nlp_solver_tol_ineq = 1e-8;
nlp_solver_tol_comp = 1e-8;
nlp_solver_ext_qp_res = 1;
qp_solver = 'partial_condensing_hpipm';
%qp_solver = 'full_condensing_hpipm';
%qp_solver = 'full_condensing_qpoases';
qp_solver_cond_N = 5;
qp_solver_cond_ric_alg = 0;
qp_solver_ric_alg = 0;
qp_solver_warm_start = 0;
qp_solver_max_iter = 100;
%sim_method = 'erk';
sim_method = 'irk';
%sim_method = 'irk_gnsf';
sim_method_num_stages = 4;
sim_method_num_steps = 3;
cost_type = 'linear_ls';
%cost_type = 'ext_cost';
model_name = 'ocp_pendulum';


%% create model entries
model = get_pendulum_on_cart_model();

% dims
nx = length(model.x);
nu = length(model.u);
ny = nu+nx; % number of outputs in lagrange term
ny_e = nx; % number of outputs in mayer term
if 0
    nbx = 0;
    nbu = nu;
    ng = 0;
    ng_e = 0;
    nh = 0;
    nh_e = 0;
else
    nbx = 0;
    nbu = 0;
    ng = 0;
    ng_e = 0;
    nh = nu;
    nh_e = 0;
end

% cost
% input-to-output matrix in lagrange term
Vu = zeros(ny, nu);
Vu(1:nu,:) = eye(nu);
% state-to-output matrix in lagrange term
Vx = zeros(ny, nx);
Vx(nu+1:end, :) = eye(nx);
% W = diag([1e-2, 1e3, 1e3, 1e-2, 1e-2]);
% high penalty on u -> no active constraints
W = diag([1e0, 1e3, 1e3, 1e-2, 1e-2]);

% terminal cost term
ny_e = nx; % number of outputs in terminal cost term
Vx_e = eye(ny_e, nx);
W_e = W(nu+1:nu+nx, nu+1:nu+nx); % weight matrix in mayer term
y_ref = zeros(ny, 1); % output reference in lagrange term
y_ref_e = zeros(ny_e, 1); % output reference in mayer term

% constraints
x0 = [0; pi; 0; 0];
%Jbx = zeros(nbx, nx); for ii=1:nbx Jbx(ii,ii)=1.0; end
%lbx = -4*ones(nbx, 1);
%ubx =  4*ones(nbx, 1);
Jbu = eye(nbu, nu);
lbu = -80*ones(nu, 1);
ubu =  80*ones(nu, 1);


%% OCP formulation
ocp = AcadosOcp();
ocp.model = model;
ocp.model.name = model_name;

if strcmp(cost_type, 'ext_cost')
    ocp.cost.cost_type = 'EXTERNAL';
    ocp.cost.cost_type_e = 'EXTERNAL';
    ocp.model.cost_expr_ext_cost = 0.5 * model.x' * diag([1e3, 1e3, 1e-2, 1e-2]) * model.x + 0.5 * model.u' * 1e-2 * model.u;
    ocp.model.cost_expr_ext_cost_e = 0.5 * model.x' * diag([1e3, 1e3, 1e-2, 1e-2]) * model.x;
else
    ocp.cost.cost_type_0 = 'LINEAR_LS';
    ocp.cost.cost_type = 'LINEAR_LS';
    ocp.cost.cost_type_e = 'LINEAR_LS';
    ocp.cost.Vu_0 = Vu;
    ocp.cost.Vx_0 = Vx;
    ocp.cost.W_0 = W;
    ocp.cost.yref_0 = y_ref;
    ocp.cost.Vu = Vu;
    ocp.cost.Vx = Vx;
    ocp.cost.Vx_e = Vx_e;
    ocp.cost.W = W;
    ocp.cost.W_e = W_e;
    ocp.cost.yref = y_ref;
    ocp.cost.yref_e = y_ref_e;
end

if strcmp(sim_method, 'erk')
    ocp.solver_options.integrator_type = 'ERK';
elseif strcmp(sim_method, 'irk_gnsf')
    ocp.solver_options.integrator_type = 'GNSF';
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
ocp.solver_options.nlp_solver_max_iter = nlp_solver_max_iter;
ocp.solver_options.nlp_solver_tol_stat = nlp_solver_tol_stat;
ocp.solver_options.nlp_solver_tol_eq = nlp_solver_tol_eq;
ocp.solver_options.nlp_solver_tol_ineq = nlp_solver_tol_ineq;
ocp.solver_options.nlp_solver_tol_comp = nlp_solver_tol_comp;
ocp.solver_options.qp_solver = upper(qp_solver);
ocp.solver_options.qp_solver_cond_N = qp_solver_cond_N;
ocp.solver_options.qp_solver_ric_alg = qp_solver_ric_alg;
ocp.solver_options.qp_solver_cond_ric_alg = qp_solver_cond_ric_alg;
ocp.solver_options.qp_solver_warm_start = qp_solver_warm_start;
ocp.solver_options.qp_solver_iter_max = qp_solver_max_iter;
ocp.solver_options.sim_method_num_stages = sim_method_num_stages;
ocp.solver_options.sim_method_num_steps = sim_method_num_steps;
ocp_solver = AcadosOcpSolver(ocp);

% set trajectory initialization
%x_traj_init = zeros(nx, N+1);
%for ii=1:N x_traj_init(:,ii) = [0; pi; 0; 0]; end
x_traj_init = [linspace(0, 0, N+1); linspace(pi, 0, N+1); linspace(0, 0, N+1); linspace(0, 0, N+1)];

u_traj_init = zeros(nu, N);

% if not set, the trajectory is initialized with the previous solution
ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);

% solve ocp
tic;
ocp_solver.solve();
time_ext = toc;

% get solution
u = ocp_solver.get('u');
x = ocp_solver.get('x');

%% evaluation
status = ocp_solver.get('status');
sqp_iter = ocp_solver.get('sqp_iter');
time_tot = ocp_solver.get('time_tot');
time_lin = ocp_solver.get('time_lin');
time_reg = ocp_solver.get('time_reg');
time_qp_sol = ocp_solver.get('time_qp_sol');

fprintf('\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n',...
    status, sqp_iter, time_ext*1e3, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3);

ocp_solver.print('stat');


%% figures
% plot trajectories
if 1
    for ii=1:N+1
        x_cur = x(:,ii);
    %    visualize;
    end

    figure;
    subplot(2,1,1);
    plot(0:N, x);
    title('trajectories')
    xlim([0 N]);
    legend('p', 'theta', 'v', 'omega');
    subplot(2,1,2);
    plot(0:N-1, u);
    xlim([0 N]);
    legend('F');
end

% plot residuals over iteraions
stat = ocp_solver.get('stat');
if 0 && (strcmp(nlp_solver, 'sqp'))
    figure;
     plot([0: size(stat,1)-1], log10(stat(:,2)), 'r-x');
     hold on
     plot([0: size(stat,1)-1], log10(stat(:,3)), 'b-x');
     plot([0: size(stat,1)-1], log10(stat(:,4)), 'g-x');
     plot([0: size(stat,1)-1], log10(stat(:,5)), 'k-x');
%    semilogy(0: size(stat,1)-1, stat(:,2), 'r-x');
%    hold on
%    semilogy(0: size(stat,1)-1, stat(:,3), 'b-x');
%    semilogy(0: size(stat,1)-1, stat(:,4), 'g-x');
%    semilogy(0: size(stat,1)-1, stat(:,5), 'k-x');
    hold off
    xlabel('iter')
    ylabel('res')
    legend('res stat', 'res eq', 'res ineq', 'res compl');
end


if status==0
    fprintf('\nsuccess!\n\n');
else
    fprintf('\nsolution failed!\n\n');
end


%% paramteric sensitivity of solution
if 1
    field = 'ex'; % equality constraint on states
    stage = 0;
    index = 1;
    ocp_solver.eval_param_sens(field, stage, index);

    sens_u = ocp_solver.get('sens_u');
    sens_x = ocp_solver.get('sens_x');

    % plot sensitivity
    figure
    subplot(2,1,1);
    plot(0:N, sens_x);
    title('sensitivities')
    xlim([0 N]);
    legend('p', 'theta', 'v', 'omega');
    subplot(2,1,2);
    plot(0:N-1, sens_u);
    xlim([0 N]);
    legend('F');

    % plot predicted solution
    figure
    subplot(2,1,1);
    plot(0:N, x+sens_x);
    title('predicted trajectories')
    xlim([0 N]);
    legend('p', 'theta', 'v', 'omega');
    subplot(2,1,2);
    plot(0:N-1, u+sens_u);
    xlim([0 N]);
    legend('F');

    for ii=1:N+1
        x_cur = x(:,ii)+sens_x(:,ii);
    %    visualize;
    end

end

sens_u = zeros(nx, N);
% get sensitivities w.r.t. initial state value with index
for index = 0:nx-1
    ocp_solver.eval_param_sens(field, stage, index);
    sens_u(index+1,:) = ocp_solver.get('sens_u');
end
disp('solution sensitivity dU_dx0')
disp(sens_u)


if is_octave()
    waitforbuttonpress;
end
