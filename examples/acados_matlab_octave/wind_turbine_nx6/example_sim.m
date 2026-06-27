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
% Wind turbine standalone simulation example ported to the NEW acados
% MATLAB/Octave interface (>= v0.4.0).

clear all
addpath('.');

% Check if env.sh is sourced
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
    disp('WARNING: env.sh does not seem to be sourced. If the build fails, run:');
    disp('source env.sh');
end

% Load simulation data
load testSim.mat

%% Arguments
method = 'IRK';
sens_forw = true;
num_stages = 4;
num_steps = 4;

%% AcadosModel
model = sim_model_wind_turbine_nx6();

nx = length(model.x);
nu = length(model.u);
np = length(model.p);

Ts = 0.2;

%% Set up AcadosSim
sim = AcadosSim();
sim.model = model;

sim.solver_options.integrator_type = method;
sim.solver_options.Tsim = Ts;
sim.solver_options.num_stages = num_stages;
sim.solver_options.num_steps = num_steps;
sim.solver_options.sens_forw = sens_forw;

% Parameter default
sim.parameter_values = zeros(np, 1);

%% Create sim solver
sim_solver = AcadosSimSolver(sim);

%% Simulation loop
% PI-controller for rotor speed tracking
uctrl = 0.0;
uctrlI = 0.0;
kI = 1e-1;
kP = 10;

nsim = 15;

x_sim = zeros(nx, nsim + 1);
x_sim(:, 1) = statesFAST(1, :);

tic;
for nn = 1:nsim

    % Compute input
    u = Usim(nn, 1:2);
    u(2) = max(u(2) - uctrl, 0);

    % Update state, input, parameter
    sim_solver.set('x', x_sim(:, nn));
    sim_solver.set('u', u);
    sim_solver.set('p', Usim(nn, 3));

    % Initial guess for implicit integrator
    if strcmp(method, 'IRK')
        sim_solver.set('xdot', zeros(nx, 1));
    end

    % Solve
    sim_solver.solve();

    x_sim(:, nn + 1) = sim_solver.get('xn');

    % Update PI controller
    ctrlErr = statesFAST(nn + 1, 1) - x_sim(1, nn + 1);
    uctrlI = uctrlI + kI * ctrlErr * Ts;
    uctrl = kP * ctrlErr + uctrlI;

end

time_solve = toc / nsim

x_sim(:, 1:nsim + 1)

% Forward sensitivities
% S_forw = sim_solver.get('S_forw');
% Sx = sim_solver.get('Sx');
% Su = sim_solver.get('Su');

fprintf('\nsuccess!\n\n');

