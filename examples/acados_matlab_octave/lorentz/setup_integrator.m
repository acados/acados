%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

function [sim_solver] = setup_integrator(model, h)
    sim = AcadosSim();
    sim.model = model;
    sim.model.name = 'lorentz_model_integrator';
    sim.solver_options.Tsim = h;

    % options
    sim.solver_options.num_stages = 2;
    sim.solver_options.num_steps = 5;
    sim.solver_options.integrator_type = 'ERK';
    sim.solver_options.sens_forw = true; % generate forward sensitivities

    % create acados integrator
    sim_solver = AcadosSimSolver(sim);
end
