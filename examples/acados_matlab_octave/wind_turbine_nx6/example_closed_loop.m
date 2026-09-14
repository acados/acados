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

%
% Wind turbine CLOSED-LOOP example (OCP + integrator) ported to the NEW
% acados MATLAB/Octave interface (>= v0.4.0).

clear all
addpath('.');

% Check if env.sh is sourced
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
    disp('WARNING: env.sh does not seem to be sourced. If the build fails, run:');
    disp('source env.sh');
end

% Get references (provides x0_ref, u0_ref, wind0_ref, y_ref)
compute_setup;

%% Arguments
% Simulation integrator
sim_method = 'IRK';
sim_sens_forw = false;
sim_num_stages = 4;
sim_num_steps = 1;

% OCP parameters
ocp_N = 40;
ocp_nlp_solver = 'SQP';
ocp_nlp_solver_exact_hessian = 'GAUSS_NEWTON';
regularize_method = 'NO_REGULARIZE';
ocp_nlp_solver_max_iter = 50;
ocp_nlp_solver_ext_qp_res = 1;
ocp_qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp_qp_solver_cond_N = 5;
ocp_qp_solver_cond_ric_alg = 0;
ocp_qp_solver_ric_alg = 0;
ocp_qp_solver_warm_start = 0;
ocp_sim_method = 'IRK';
ocp_sim_method_num_stages = 4;
ocp_sim_method_num_steps = 1;
cost_type = 'LINEAR_LS';

%% Create model entries
model = ocp_model_wind_turbine_nx6();

nx = length(model.x);   % 8
nu = length(model.u);   % 2
np = length(model.p);   % 1

%% Dimensions
Ts = 0.2; % Sampling time
T = ocp_N * Ts; % Horizon length [s]
ny = 4; % Outputs in Lagrange term
ny_e = 2; % Outputs in Mayer term

% Soft-constraint counts
ns = 2;
ns_e = 2;

%% Set up OCP
ocp = AcadosOcp();
ocp.model = model;

%% Cost (Linear least squares)
% State-to-output matrix (Lagrange)
Vx = zeros(ny, nx);
Vx(1, 1) = 1.0;
Vx(2, 5) = 1.0;

% Input-to-output matrix (Lagrange)
Vu = zeros(ny, nu);
Vu(3, 1) = 1.0;
Vu(4, 2) = 1.0;

% State-to-output matrix (Mayer)
Vx_e = zeros(ny_e, nx);
Vx_e(1, 1) = 1.0;
Vx_e(2, 5) = 1.0;

% Weight matrix (Lagrange)
W = zeros(ny, ny);
W(1, 1) = 1.5114;
W(2, 1) = -0.0649;
W(1, 2) = -0.0649;
W(2, 2) = 0.0180;
W(3, 3) = 0.01;
W(4, 4) = 0.001;

% Weight matrix (Mayer)
W_e = zeros(ny_e, ny_e);
W_e(1, 1) = 1.5114;
W_e(2, 1) = -0.0649;
W_e(1, 2) = -0.0649;
W_e(2, 2) = 0.0180;

ocp.cost.cost_type = cost_type;
ocp.cost.cost_type_e = cost_type;
ocp.cost.cost_type_0 = cost_type;

ocp.cost.Vx_0 = Vx;
ocp.cost.Vu_0 = Vu;
ocp.cost.W_0 = W;
ocp.cost.yref_0 = zeros(ny, 1);

ocp.cost.Vx = Vx;
ocp.cost.Vu = Vu;
ocp.cost.W = W;
ocp.cost.yref = zeros(ny, 1);

ocp.cost.Vx_e = Vx_e;
ocp.cost.W_e = W_e;
ocp.cost.yref_e = zeros(ny_e, 1);

% Slack penalties
ocp.cost.Zl = 1e2 * ones(ns, 1);
ocp.cost.Zu = 1e2 * ones(ns, 1);
ocp.cost.zl = zeros(ns, 1);
ocp.cost.zu = zeros(ns, 1);

ocp.cost.Zl_e = 1e2 * ones(ns_e, 1);
ocp.cost.Zu_e = 1e2 * ones(ns_e, 1);
ocp.cost.zl_e = zeros(ns_e, 1);
ocp.cost.zu_e = zeros(ns_e, 1);

%% Constraints
% Bound constants
dbeta_min = -8.0;
dbeta_max = 8.0;
dM_gen_min = -1.0;
dM_gen_max = 1.0;
OmegaR_min = 6.0 / 60 * 2 * pi;
OmegaR_max = 13.0 / 60 * 2 * pi;
beta_min = 0.0;
beta_max = 35.0;
M_gen_min = 0.0;
M_gen_max = 5.0;
Pel_min = 0.0;
Pel_max = 5.0;

% State bounds
ocp.constraints.idxbx = [0; 6; 7];
ocp.constraints.lbx = [OmegaR_min; beta_min; M_gen_min];
ocp.constraints.ubx = [OmegaR_max; beta_max; M_gen_max];

% Terminal state bounds
ocp.constraints.idxbx_e = [0; 6; 7];
ocp.constraints.lbx_e = [OmegaR_min; beta_min; M_gen_min];
ocp.constraints.ubx_e = [OmegaR_max; beta_max; M_gen_max];

% Input bounds
ocp.constraints.idxbu = [0; 1];
ocp.constraints.lbu = [dbeta_min; dM_gen_min];
ocp.constraints.ubu = [dbeta_max; dM_gen_max];

% Nonlinear power constraint
ocp.constraints.lh = Pel_min;
ocp.constraints.uh = Pel_max;
ocp.constraints.lh_e = Pel_min;
ocp.constraints.uh_e = Pel_max;

% Soft box state constraint
ocp.constraints.idxsbx = 0;
ocp.constraints.idxsbx_e = 0;

% Soft nonlinear constraint
ocp.constraints.idxsh = 0;
ocp.constraints.idxsh_e = 0;

% Initial state
ocp.constraints.x0 = x0_ref;

% Parameter default
ocp.parameter_values = zeros(np, 1);

%% Create OCP solver
ocp.solver_options.N_horizon = ocp_N;
ocp.solver_options.tf = T;

ocp.solver_options.nlp_solver_type = ocp_nlp_solver;
ocp.solver_options.hessian_approx = ocp_nlp_solver_exact_hessian;
ocp.solver_options.regularize_method = regularize_method;
ocp.solver_options.nlp_solver_ext_qp_res = ocp_nlp_solver_ext_qp_res;
ocp.solver_options.nlp_solver_max_iter = ocp_nlp_solver_max_iter;

ocp.solver_options.qp_solver = ocp_qp_solver;
ocp.solver_options.qp_solver_cond_N = ocp_qp_solver_cond_N;
ocp.solver_options.qp_solver_cond_ric_alg = ocp_qp_solver_cond_ric_alg;
ocp.solver_options.qp_solver_ric_alg = ocp_qp_solver_ric_alg;
ocp.solver_options.qp_solver_warm_start = ocp_qp_solver_warm_start;

ocp.solver_options.integrator_type = ocp_sim_method;
ocp.solver_options.sim_method_num_stages = ocp_sim_method_num_stages;
ocp.solver_options.sim_method_num_steps = ocp_sim_method_num_steps;

ocp_solver = AcadosOcpSolver(ocp);

%% AcadosSim for the plant integrator
sim = AcadosSim();
sim.model = model;
sim.solver_options.integrator_type = sim_method;
sim.solver_options.Tsim = T / ocp_N;
sim.solver_options.num_stages = sim_num_stages;
sim.solver_options.num_steps = sim_num_steps;
sim.solver_options.sens_forw = sim_sens_forw;
sim.parameter_values = zeros(np, 1);

sim_solver = AcadosSimSolver(sim);

%% Closed-loop simulation
n_sim = 100;
n_sim_max = length(wind0_ref) - ocp_N;
if n_sim > n_sim_max
    n_sim = n_sim_max;
end

x_sim = zeros(nx, n_sim + 1);
x_sim(:, 1) = x0_ref; % Initial state
u_sim = zeros(nu, n_sim);

sqp_iter_sim = zeros(n_sim, 1);

% Trajectory initialization
x_traj_init = repmat(x0_ref, 1, ocp_N + 1);
u_traj_init = repmat(u0_ref, 1, ocp_N);
pi_traj_init = zeros(nx, ocp_N);

for ii = 1:n_sim
    tic

    % Set x0
    ocp_solver.set('constr_x0', x_sim(:, ii));

    % Set parameter (wind) per stage
    for jj = 0:ocp_N - 1
        ocp_solver.set('p', wind0_ref(:, ii + jj), jj);
    end

    % Set reference per stage
    for jj = 0:ocp_N - 1
        ocp_solver.set('cost_y_ref', y_ref(:, ii + jj), jj);
    end
    ocp_solver.set('cost_y_ref', y_ref(1:ny_e, ii + ocp_N), ocp_N);

    % Trajectory initialization
    ocp_solver.set('init_x', x_traj_init);
    ocp_solver.set('init_u', u_traj_init);
    ocp_solver.set('init_pi', pi_traj_init);

    % Solve
    ocp_solver.solve();

    % Get solution
    x = ocp_solver.get('x');
    u = ocp_solver.get('u');
    pi = ocp_solver.get('pi');

    % Store first input
    u_sim(:, ii) = ocp_solver.get('u', 0);

    % Integrate plant one step
    sim_solver.set('x', x_sim(:, ii));
    sim_solver.set('u', u_sim(:, ii));
    sim_solver.set('p', wind0_ref(:, ii));
    sim_solver.solve();
    x_sim(:, ii + 1) = sim_solver.get('xn');

    % Shift trajectory for next initialization
    x_traj_init = [x(:, 2:ocp_N + 1), x(:, ocp_N + 1)];
    u_traj_init = [u(:, 2:ocp_N), u(:, ocp_N)];
    pi_traj_init = [pi(:, 2:ocp_N), pi(:, ocp_N)];

    time_ext = toc;

    electrical_power = 0.944 * 97 / 100 * x(1, 1) * x(6, 1);

    status = ocp_solver.get('status');
    sqp_iter = ocp_solver.get('sqp_iter');
    time_tot = ocp_solver.get('time_tot');
    time_lin = ocp_solver.get('time_lin');
    time_qp_sol = ocp_solver.get('time_qp_sol');

    sqp_iter_sim(ii) = sqp_iter;

    fprintf(['\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] ', ...
             '(time_lin = %f [ms], time_qp_sol = %f [ms]), Pel = %f\n'], ...
            status, sqp_iter, time_ext * 1e3, time_tot * 1e3, ...
            time_lin * 1e3, time_qp_sol * 1e3, electrical_power);

    if 0
        ocp_solver.print('stat')
    end
end

electrical_power = 0.944 * 97 / 100 * x_sim(1, :) .* x_sim(6, :);

if status == 0
    fprintf('\nsuccess!\n\n');
else
    fprintf('\nsolution failed!\n\n');
end

%% Figures
if 1
    figure;
    
    subplot(3, 1, 1);
    plot(0:n_sim, x_sim);
    xlim([0, n_sim]);
    ylabel('states');
    
    subplot(3, 1, 2);
    plot(0:n_sim - 1, u_sim);
    xlim([0, n_sim]);
    ylabel('inputs');
    
    subplot(3, 1, 3);
    plot(0:n_sim, electrical_power);
    hold on;
    plot([0, n_sim], [Pel_max, Pel_max]);
    hold off;
    xlim([0, n_sim]);
    ylim([4.0, 6.0]);
    ylabel('electrical power');

    figure;
    plot(1:n_sim, sqp_iter_sim, 'rx');
    hold on;
    plot([1, n_sim], [ocp_nlp_solver_max_iter, ocp_nlp_solver_max_iter]);
    hold off;
    ylim([0, ocp_nlp_solver_max_iter + 1]);
    ylabel('sqp iterations');
    xlabel('sqp calls');
    
    if is_octave()
        waitforbuttonpress;
    end
end
