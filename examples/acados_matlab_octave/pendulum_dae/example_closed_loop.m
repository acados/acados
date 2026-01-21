%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;

clear all; clc;

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
    error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end

%% model
model = get_pendulum_dae_model();

nx = length(model.x);
nu = length(model.u);
nz = length(model.z);
ny = nu+nx; % number of outputs in lagrange term
ny_e = nx; % number of outputs in mayer term
ng = 0; % number of general linear constraints intermediate stages
nbx = 0; % number of bounds on state x
nbu = nu; % number of bounds on controls u

%% Acados Ocp
ocp = AcadosOcp();

ocp_N = 50;
h = 0.05;
T = ocp_N*h;

ocp.model = model;
ocp.solver_options.tf = T;
ocp.solver_options.N_horizon = ocp_N;
ocp.solver_options.nlp_solver_type = 'SQP_RTI'; % 'SQP', 'SQP_RTI'
ocp.solver_options.hessian_approx = 'EXACT'; % 'EXACT', 'GAUSS_NEWTON'
ocp.solver_options.regularize_method = 'PROJECT_REDUC_HESS';% NO_REGULARIZE, PROJECT, PROJECT_REDUC_HESS, MIRROR, CONVEXIFY
ocp.solver_options.nlp_solver_max_iter = 100;
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_cond_N = 5; % for partial condensing
ocp.solver_options.integrator_type = 'IRK'; % 'DISCRETE','ERK','IRK
ocp.solver_options.qp_solver_cond_ric_alg = 0;
ocp.solver_options.qp_solver_ric_alg = 0;
ocp.solver_options.qp_solver_warm_start = 0;
ocp.solver_options.sim_method_num_stages = 4 * ones(ocp_N, 1);
ocp.solver_options.sim_method_num_steps = 1;
ocp.solver_options.sim_method_newton_iter = 3;
ocp.solver_options.compile_interface = [];

%  sym_x = [xpos, ypos, alpha, vx, vy, valpha]
length_pendulum = 5;
xsteady = [ 0; -length_pendulum; 0; 0; 0; 0];
xtarget = [ 0; +length_pendulum; pi; 0; 0; 0];
% xtarget = xsteady;
uref = 0;

alpha0 = -.01;
xp0 = length_pendulum * sin(alpha0);
yp0 = - length_pendulum * cos(alpha0);
x0 = [ xp0; yp0; alpha0; 0; 0; 0];
% x0 = xsteady + 1e-4 * ones(nx,1);

disp('state')
disp(model.x)

% cost
% linear least square cost: y^T * W * y, where y = Vx * x + Vu * u - y_ref
Vx = eye(ny, nx); % state-to-output matrix in lagrange term
Vu = zeros(ny, nu);
Vu(nx:end, nx:end) = eye(nu); % input-to-output matrix in lagrange term
Vx_e = Vx(1:ny_e,:); % state-to-output matrix in mayer term
% weight matrix in lagrange term
W = diag([1e3 * ones(3,1);... % xpos, ypos, alpha
            ones(3,1); ... %speeds
            1e2 ]); % control
W_e = W(nu+1:nu+nx, nu+1:nu+nx); % weight matrix in mayer term
yr = [xtarget; uref]; % output reference in lagrange term
yr_e = xtarget; % output reference in mayer term

ocp.cost.cost_type = 'LINEAR_LS';
ocp.cost.cost_type_e = 'LINEAR_LS';

ocp.cost.Vu = Vu;
ocp.cost.Vx = Vx;
ocp.cost.Vx_e = Vx_e;
ocp.cost.W = W;
ocp.cost.W_e = W_e;
ocp.cost.yref = yr;
ocp.cost.yref_e = yr_e;
ocp.cost.Vz = zeros(ny,nz);

% constraints
ocp.constraints.x0 = x0;
constraint_h = 1;
if constraint_h
    nh = length(model.con_h_expr);
    ocp.model.con_h_expr_0 = model.con_h_expr;
    ocp.model.con_h_expr = model.con_h_expr;
    lh =    0;
    uh =  200;
    ocp.constraints.lh_0 = lh;
    ocp.constraints.uh_0 = uh;
    ocp.constraints.lh = lh;
    ocp.constraints.uh = uh;
else
    nh = 0;
end

lbu = -80*ones(nu, 1);
ubu =  80*ones(nu, 1);
ocp.constraints.idxbu = (0:nbu-1)';
ocp.constraints.lbu = lbu;
ocp.constraints.ubu = ubu;
%% acados ocp
ocp_solver = AcadosOcpSolver(ocp);

%% acados sim model
sim = AcadosSim();
sim.model = model;
sim.solver_options.Tsim = h;
sim.solver_options.integrator_type = 'IRK';  % 'ERK', 'IRK'
sim.solver_options.sens_forw = false; % true, false
sim.solver_options.jac_reuse = false; % true, false
sim.solver_options.num_stages = 3;
sim.solver_options.num_steps = 3;
sim.solver_options.newton_iter = 3;
sim.solver_options.compile_interface = 'AUTO';

%% acados sim solver
% create integrator
sim_solver = AcadosSimSolver(sim);

%% closed loop simulation
N_sim = 99;

x_sim = zeros(nx, N_sim+1);
x_sim(:,1) = x0; % initial state
u_sim = zeros(nu, N_sim);
z_sim = zeros(nz, N_sim);

% initialization
xdot0 = zeros(nx, 1);
z0 = zeros(nz, 1);

% set trajectory initialization
x_traj_init = repmat(x0, 1, ocp_N + 1);
u_traj_init = zeros(nu, ocp_N);
pi_traj_init = zeros(nx, ocp_N);
z_traj_init = 0*ones(nz, ocp_N);
xdot_traj_init = 0*ones(nx, ocp_N);

sqp_iter = zeros(N_sim,1);
sqp_time = zeros(N_sim,1);

tic
for ii=1:N_sim

    % set initial state
    ocp_solver.set('constr_x0', x_sim(:,ii));
    % set trajectory initialization (if not, set internally using previous solution)
    ocp_solver.set('init_x', x_traj_init);
    ocp_solver.set('init_u', u_traj_init);
    ocp_solver.set('init_pi', pi_traj_init);
    if ii == 1
        ocp_solver.set('init_z', z_traj_init);
        ocp_solver.set('init_xdot', xdot_traj_init);
    end

    ocp_solver.solve();

    status = ocp_solver.get('status');
    sqp_iter(ii) = ocp_solver.get('sqp_iter');
    sqp_time(ii) = ocp_solver.get('time_tot');
    if status ~= 0
        ocp.print('stat');
        error('solver returned nonzero exit status.')
    end

    % get solution for initialization of next NLP
    x_traj = ocp_solver.get('x');
    u_traj = ocp_solver.get('u');
    pi_traj = ocp_solver.get('pi');
    z_traj = ocp_solver.get('z');

    % shift trajectory for initialization
    x_traj_init = [x_traj(:,2:end), x_traj(:,end)];
    u_traj_init = [u_traj(:,2:end), u_traj(:,end)];
    pi_traj_init = [pi_traj(:,2:end), pi_traj(:,end)];
    z_traj_init = [z_traj(:,2:end), z_traj(:,end)];

    u_sim(:,ii) = ocp_solver.get('u', 0); % get control input
    sim_solver.set('xdot', xdot0);
    sim_solver.set('z', z0);
    sim_solver.set('x', x_sim(:,ii));     % set initial state
    sim_solver.set('u', u_sim(:,ii));     % set input
    sim_solver.solve();    % simulate state

    % get simulated state
    x_sim(:,ii+1) = sim_solver.get('xn');
    z_sim(:,ii) = sim_solver.get('zn');

end
format short e

toc

% trajectory plot
figure;
subplot(4, 1, 1);
plot(1:N_sim+1, x_sim(1:2,:));
legend('xpos', 'ypos');

subplot(4, 1, 2);
plot(1:N_sim+1, x_sim(3,:));
legend('alpha')

subplot(4, 1, 3);
plot(1:N_sim+1, x_sim(4:6,:));
legend('vx', 'vy', 'valpha');

subplot(4,1,4);
plot(1:N_sim+1, [0 u_sim]);
legend('u')

if is_octave()
    waitforbuttonpress;
end


% iterations, CPU time
figure();
subplot(2, 1, 1);
plot(1:N_sim, sqp_iter)
ylabel('# iterations')
xlabel('ocp instance')
subplot(2, 1, 2);
plot(1:N_sim, sqp_time)
ylabel('CPU time')
xlabel('ocp instance')
if is_octave()
    waitforbuttonpress;
end

% check consistency
xp = x_sim(1,:);
yp = x_sim(2,:);
check = abs(xp.^2 + yp.^2 - length_pendulum^2);
tol_pendulum = 1e-4;

disp(['checking for constant pendulum length, got ' num2str(check)])
if any( max(abs(check)) > tol_pendulum )
    error(['note: check for constant pendulum length failed, violation >' ...
        num2str(tol_pendulum)]);
end

% eval constraint h
ax_ = z_sim(1,:);
ay_ = z_sim(2,:);
h_vals = ax_.^2 + ay_.^2;
h_violations = [ uh - h_vals; h_vals - lh ];
max_h_violation = max(abs( h_violations( h_violations < 0 )));
if isempty(max_h_violation)
    max_h_violation = 0;
end
disp(['maximal constraint h violation   ' num2str( max_h_violation, '%e' ) ])

