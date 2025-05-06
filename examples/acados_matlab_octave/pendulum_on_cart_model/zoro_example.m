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
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS'
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

import casadi.*


%
check_acados_requirements()

%% solver settings
N = 20; % number of discretization steps
T = 1; % [s] prediction horizon length
x0 = [0.0, 0.15*pi, 0.0, 0.0];


%% model dynamics
model = get_pendulum_on_cart_model();
nx = length(model.x); % state size
nu = length(model.u); % input size

%% OCP formulation object
ocp = AcadosOcp();
ocp.model = model;

%% cost in nonlinear least squares form
W_x = diag([1e3, 1e3, 1e-2, 1e-2]);
W_u = 1e-2;

% initial cost term
ocp.cost.cost_type_0 = 'NONLINEAR_LS';
ocp.cost.W_0 = W_u;
ocp.cost.yref_0 = zeros(nu, 1);
ocp.model.cost_y_expr_0 = model.u;

% path cost term
ocp.cost.cost_type = 'NONLINEAR_LS';
ocp.cost.W = blkdiag(W_x, W_u);
ocp.cost.yref = zeros(nx+nu, 1);
ocp.model.cost_y_expr = vertcat(model.x, model.u);

% terminal cost term
ocp.cost.cost_type_e = 'NONLINEAR_LS';
ocp.model.cost_y_expr_e = model.x;
ocp.cost.yref_e = zeros(nx, 1);
ocp.cost.W_e = W_x;

%% define constraints
% bound on u
Fmax = 40;
ocp.constraints.lbu = [-Fmax];
ocp.constraints.ubu = [+Fmax];
ocp.constraints.idxbu = [0];

% bound on x
theta_min = -pi * 0.15;
theta_max = pi * 0.3;
ocp.constraints.lbx = [theta_min];
ocp.constraints.ubx = [theta_max];
ocp.constraints.idxbx = [1];

% bound on the terminal state
ocp.constraints.lbx_e = [theta_min];
ocp.constraints.ubx_e = [theta_max];
ocp.constraints.idxbx_e = [1];

% initial state
ocp.constraints.x0 = x0;

% define solver options
ocp.solver_options.N_horizon = N;
ocp.solver_options.tf = 1.0;
ocp.solver_options.nlp_solver_type = 'SQP_RTI';
ocp.solver_options.integrator_type = 'ERK';
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';


% custom update: disturbance propagation
ocp.solver_options.custom_update_filename = 'custom_update_function.c';
ocp.solver_options.custom_update_header_filename = 'custom_update_function.h';

ocp.solver_options.custom_update_copy = true;
ocp.solver_options.custom_templates = {
    {'custom_update_function_zoro_template.in.c', 'custom_update_function.c'},
    {'custom_update_function_zoro_template.in.h', 'custom_update_function.h'},
};

% zoro stuff
zoro_description = ZoroDescription();
zoro_description.backoff_scaling_gamma = 2;
zoro_description.P0_mat = zeros(nx, nx);
zoro_description.fdbk_K_mat = [[0.0, 0.0, 10.0, 10.0]];
zoro_description.W_mat = diag([5*1e-6, 5*1e-6, 1*1e-4, 1*1e-4]);
zoro_description.idx_lbu_t = [0];
zoro_description.idx_ubu_t = [0];
zoro_description.idx_lbx_t = [0];
zoro_description.idx_lbx_e_t = [0];
ocp.zoro_description = zoro_description;

% create solver
ocp_solver = AcadosOcpSolver(ocp);

Nsim = 100;
simX = zeros(Nsim+1, nx);
simU = zeros(Nsim, nu);
simX(1,:) = x0;

% zoro parameters
max_zoro_iter = 100;
zoro_tol = 1e-6;

% sample disturbances
% NOTE: this requires statistics addon in MATLAB and octave -> just load
% some data instead
% dist_samples = mvnrnd(zeros(nx, 1), zoro_description.W_mat, Nsim);
load('dist_samples.mat')

for idx_sim = 1:Nsim
    residuals = inf;
    % set initial state
    ocp_solver.set('constr_lbx', simX(idx_sim,:), 0)
    ocp_solver.set('constr_ubx', simX(idx_sim,:), 0)

    for idx_iter = 1:max_zoro_iter
        % preparation phase
        ocp_solver.set('rti_phase', 1)
        ocp_solver.solve();
        % constraint tightening
        ocp_solver.custom_update([]);
        % call SQP_RTI solver: feedback phase
        ocp_solver.set('rti_phase', 2)
        ocp_solver.solve()

        status = ocp_solver.get('status');

        if status ~= 0
            disp(['simU(idx_sim, :) = ', mat2str(simU(idx_sim, :))]);
            error('acados returned status %d at idx_sim = %d for initial state %s.', status, idx_sim, mat2str(simX(idx_sim, :)));
        end

        residuals = ocp_solver.get('residuals');
        if max(residuals) <= zoro_tol
            break
        end
    end

    if max(residuals) > zoro_tol
        disp('zoro does not converge');
    end

    % get solution
    simU(idx_sim,:) = ocp_solver.get('u', 0);
    simX(idx_sim+1,:) = ocp_solver.get('x', 1) + dist_samples(idx_sim);

    if idx_sim == 1
        % print backoffs
        disp(['backoff in the terminal state constraint=', ...
            num2str((theta_max - theta_min) - ocp_solver.get('qp_ubx', N) - ocp_solver.get('qp_lbx', N))])
        disp(['backoff in the N-1 input constraint=', ...
            num2str(2 * Fmax - (ocp_solver.get('qp_ubu', N-1) - ocp_solver.get('qp_lbu', N-1)))])
        % print costs
        cost = ocp_solver.get_cost();
        disp(['cost function value of solution = ', num2str(cost)]);
    end
end


%% plots
h = T / N;
ts = linspace(0, Nsim*h, Nsim+1);
figure; hold on;
States = {'p', 'theta', 'v', 'dtheta'};
for i=1:length(States)
    subplot(length(States), 1, i);
    plot(ts, simX(:,i)); grid on;
    ylabel(States{i});
    xlabel('t [s]')
end

figure
stairs(ts, [simU; simU(end)])
ylabel('F [N]')
xlabel('t [s]')
grid on

if is_octave
    waitforbuttonpress;
end
