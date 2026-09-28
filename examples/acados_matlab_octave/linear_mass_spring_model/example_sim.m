%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

clear all; clc;
check_acados_requirements();

%% model
model = get_linear_mass_spring_model();
nx = length(model.x);
nu = length(model.u);
%% set up sim
sim = AcadosSim();
sim.model = model;
sim.solver_options.Tsim = 0.5;
sim.solver_options.integrator_type = 'IRK';  % 'ERK', 'IRK'
sim.solver_options.num_stages = 4;
sim.solver_options.num_steps = 4;
sim.solver_options.sens_forw = true;
%% acados sim
% create sim
sim_solver = AcadosSimSolver(sim);
% (re)set numerical part of model
%sim_solver.set('T', 0.5);

x0 = ones(nx, 1); %x0(1) = 2.0;
tic;
sim_solver.set('x', x0);
time_set_x = toc

u = ones(nu, 1);
sim_solver.set('u', u);

% solve
tic;
sim_solver.solve();
time_solve = toc


% xn
xn = sim_solver.get('xn');
xn
% S_forw
S_forw = sim_solver.get('S_forw');
S_forw
% Sx
Sx = sim_solver.get('Sx');
Sx
% Su
Su = sim_solver.get('Su');
Su


fprintf('\nsuccess!\n\n');

