%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

clear VARIABLES

addpath('../linear_mass_spring_model/');

method = 'IRK';
sens_forw = true;
num_stages = 4;
num_steps = 4;

Ts = 0.1;

model = get_linear_mass_spring_model();

model_name = ['lin_mass_' method];
nx = length(model.x);
nu = length(model.u);

sim = AcadosSim();
sim.model = model;
sim.model.name = model_name;
sim.solver_options.Tsim = Ts;
sim.solver_options.integrator_type = method;
sim.solver_options.num_stages = num_stages;
sim.solver_options.num_steps = num_steps;
sim.solver_options.sens_forw = sens_forw;
sim_solver = AcadosSimSolver(sim);

try
    sim_solver.set('x', zeros(nx+1, 1));
    error('test_checks: setter accepted a state with the wrong dimension');
catch exception
    if ~isempty(strfind(exception.message, 'wrong dimension'))
        disp('Success: setter rejects a state with the wrong dimension')
    else
        rethrow(exception);
    end
end
