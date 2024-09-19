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
% author: Katrin Baumgaertner

clear all

h = 0.05;   % time step
N_mhe = 15;     % MHE horizon

%% model
model = lorentz_model();

sim_solver = setup_integrator(model, h);
estimator = setup_estimator(model, h, N_mhe);

nx = estimator.ocp.dims.nx;
nu = estimator.ocp.dims.nu;
ny = estimator.ocp.dims.ny - estimator.ocp.dims.nu; % cost ny = dimension of measurement y + noise w

%% Simulation
N_sim = 120;

iter_step = 50;
step = 3;

x0 = [-1 3 4 25];
x0_bar = [-1 3 4 25];

v_std = 0.0;   % standard deviation of measurement noise

x_sim = zeros(nx, N_sim+1);
y_sim = zeros(ny, N_sim);

x_sim(:,1) = x0;

for n=1:N_sim

    % set initial state
    sim_solver.set('x', x_sim(:,n));

    % solve
    sim_solver.solve();

    % get simulated state
    x_sim(:,n+1) = sim_solver.get('xn');

    % unmodeled step change in x(4)
    if n == iter_step
        x_sim(end, n+1) = x_sim(end, n+1) + step;
    end

    % measurement
    y_sim(:, n) = x_sim(1, n) + v_std*randn(1, 1);
end

%% Estimation
x_est = zeros(nx, N_sim-N_mhe);

yref_0 = zeros(ny + nu + nx, 1);
yref = zeros(ny + nu, 1);

for n=1:N_sim-N_mhe

    % set measurements
    yref_0(1:ny) = y_sim(:, n);
    yref_0(ny+nu+1:end) = x0_bar;

    estimator.set('cost_y_ref', yref_0, 0);

    for i=1:N_mhe-1
        yref(1:ny) = y_sim(:, n+i);
        estimator.set('cost_y_ref', yref, i);
    end

    %estimator.set('init_x', x_sim(:, n:n+N_mhe))
    % solve
    estimator.solve()

    x_est(:, n) = estimator.get('x', N_mhe);

    % update arrival cost (TODO: update P0 as well)
    x0_bar = estimator.get('x', 1);
end

%% Plot
ts = h*(0:N_sim);

figure;
States = {'x_1', 'x_2', 'x_3', 'p'};
for i=1:length(States)
    subplot(length(States), 1, i); hold on;
    plot(ts, x_sim(i,:));
    plot(ts(N_mhe+1:end-1), x_est(i,:));

    if i == 1
        plot(ts(1:end-1), y_sim, 'x');
        legend('true', 'est', 'measured');
    end
    grid on;
    ylabel(States{i});
    xlabel('t [s]');
end

figure;
States = {'abs. error x_1', 'abs. error x_2', 'abs. error x_3', 'abs. error p'};
for i=1:length(States)
    subplot(length(States), 1, i); hold on; grid on;

    plot(ts(N_mhe+1:end-1), abs(x_est(i,:) - x_sim(i, N_mhe+1:end-1)));

    ylabel(States{i});
    xlabel('t [s]');
end
