%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

clear all

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end

%% arguments
compile_interface = 'auto';
gnsf_detect_struct = 'true';
%method = 'erk';
% method = 'irk';
method = 'irk_gnsf';
sens_forw = 'true';
jac_reuse = 'true';
num_stages = 4;
num_steps = 4;
newton_iter = 5;
model_name = 'sim_pendulum';

h = 0.1;
x0 = [0; 1e-1; 0; 0e0];
u = 0;

%% model
model = get_pendulum_on_cart_model();

nx = length(model.x);
nu = length(model.u);

%% Simulation formulation
sim = AcadosSim();
sim.model = model;
sim.model.name = model_name;
sim.solver_options.Tsim = h;
if strcmp(method, 'irk_gnsf')
    sim.solver_options.integrator_type = 'GNSF';
else
    sim.solver_options.integrator_type = upper(method);
end
sim.solver_options.num_stages = num_stages;
sim.solver_options.num_steps = num_steps;
sim.solver_options.newton_iter = newton_iter;
sim.solver_options.sens_forw = strcmp(sens_forw, 'true');
sim.solver_options.jac_reuse = strcmp(jac_reuse, 'true');
sim_solver = AcadosSimSolver(sim);
% (re)set numerical part of model
%sim_solver.set('T', 0.5);
%sim_solver.C_sim
%sim_solver.C_sim_ext_fun


N_sim = 100;

x_sim = zeros(nx, N_sim+1);
x_sim(:,1) = x0;

tic
for ii=1:N_sim

	% set initial state
	sim_solver.set('x', x_sim(:,ii));
	sim_solver.set('u', u);

    % initialize implicit integrator
    if (strcmp(method, 'irk'))
        sim_solver.set('xdot', zeros(nx,1));
    elseif (strcmp(method, 'irk_gnsf'))
        n_out = sim_solver.sim.model.gnsf_model.dims.nout;
        sim_solver.set('phi_guess', zeros(n_out,1));
    end

	% solve
	sim_solver.solve();


	% get simulated state
	x_sim(:,ii+1) = sim_solver.get('xn');

end
simulation_time = toc


% xn
xn = sim_solver.get('xn');
xn
% S_forw
S_forw = sim_solver.get('S_forw')
Sx = sim_solver.get('Sx');
Su = sim_solver.get('Su');

%x_sim

% for ii=1:N_sim+1
% 	x_cur = x_sim(:,ii);
% 	visualize;
% end

figure;
plot(1:N_sim+1, x_sim);
legend('p', 'theta', 'v', 'omega');


fprintf('\nsuccess!\n\n');


if is_octave()
    waitforbuttonpress;
end
