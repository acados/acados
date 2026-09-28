%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

%% minimal example of acados integrator matlab interface
clear all; clc;

addpath('../pendulum_on_cart_model')

check_acados_requirements()


% simulation parameters
N_sim = 100;
x0 = [0; 1e-1; 0; 0]; % initial state
u0 = 0; % control input

%% define model dynamics
model = get_pendulum_on_cart_model();
nx = length(model.x);

sim = AcadosSim();
sim.model = model;
sim.solver_options.Tsim = 0.1; % simulation time
sim.solver_options.integrator_type = 'ERK';

%% create integrator
sim_solver = AcadosSimSolver(sim);

%% simulate system in loop
x_sim = zeros(nx, N_sim+1);
x_sim(:,1) = x0;

for ii=1:N_sim

    % set initial state
    sim_solver.set('x', x_sim(:,ii));
    sim_solver.set('u', u0);

    % initialize implicit integrator
    if strcmp(sim.solver_options.integrator_type, 'IRK')
        sim_solver.set('xdot', zeros(nx,1));
    elseif strcmp(sim.solver_options.integrator_type, 'GNSF')
        n_out = sim_solver.sim.model.gnsf_model.dims.nout;
        sim_solver.set('phi_guess', zeros(n_out,1));
    end

    % solve
    sim_solver.solve();

    % get simulated state
    x_sim(:,ii+1) = sim_solver.get('xn');

    % one can also use:
    % x_sim(:,ii+1) = sim_solver.simulate(x_sim(:,ii), u0);
end

% forward sensitivities ( dxn_d[x0,u] )
S_forw = sim_solver.get('S_forw');

%% plot state trajectories
ts = linspace(0, sim.solver_options.Tsim*(N_sim+1), N_sim+1);
figure; hold on;
States = {'p', 'theta', 'v', 'dtheta'};
for i=1:length(States)
    subplot(length(States), 1, i);
    plot(ts, x_sim(i,:)); grid on;
    ylabel(States{i});
    xlabel('t [s]')
end