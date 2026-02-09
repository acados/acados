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
disp('state')
disp(model.x)

nx = length(model.x);
nu = length(model.u);
nz = length(model.z);

%% acados sim model
sim = AcadosSim();
sim.model = model;
sim.solver_options.Tsim = 0.1;
sim.solver_options.integrator_type = 'IRK';  % 'ERK', 'IRK'
sim.solver_options.sens_forw = true; % true, false
sim.solver_options.jac_reuse = false; % true, false
sim.solver_options.num_stages = 3;
sim.solver_options.num_steps = 3;
sim.solver_options.newton_iter = 3;
sim.solver_options.compile_interface = 'AUTO';

length_pendulum = 5;

alpha0 = 0.1;
xp0 = length_pendulum * sin(alpha0);
yp0 = - length_pendulum * cos(alpha0);
x0 = [ xp0; yp0; alpha0; 0; 0; 0];

u = 0;

% steady state
% x0 = [ 0; -length_pendulum; 0; 0; 0; 0];

%% acados sim
% create integrator
sim_solver = AcadosSimSolver(sim);

N_sim = 100;

x_sim = zeros(nx, N_sim+1);
x_sim(:,1) = x0;

% initialization
xdot0 = zeros(nx, 1);
z0 = zeros(nz, 1);

tic
for ii=1:N_sim

	% set initial state
	sim_solver.set('x', x_sim(:,ii));
	sim_solver.set('u', u);

	% set adjoint seed
	sim_solver.set('seed_adj', ones(nx,1));
    sim_solver.set('xdot', xdot0);
    sim_solver.set('z', z0);
	% solve
	sim_solver.solve();

	% get simulated state
	x_sim(:,ii+1) = sim_solver.get('xn');

	% forward sensitivities
	S_forw = sim_solver.get('S_forw');
	Sx = sim_solver.get('Sx');
	Su = sim_solver.get('Su');

end
format short e
xfinal = sim_solver.get('xn')';
S_adj = sim_solver.get('S_adj')';
z = sim_solver.get('zn')'; % approximate value of algebraic variables at start of simulation
S_alg = sim_solver.get('S_algebraic'); % sensitivities of algebraic variables z

simulation_time = toc


figure;
Nsub = 4;
subplot(Nsub, 1, 1);
plot(1:N_sim+1, x_sim(1,:));
legend('xpos');

subplot(Nsub, 1, 2);
plot(1:N_sim+1, x_sim(2,:));
legend('ypos');

subplot(Nsub, 1, 3);
plot(1:N_sim+1, x_sim(3,:));
legend('alpha')

subplot(Nsub, 1, 4);
plot(1:N_sim+1, x_sim(4:6,:));
legend('vx', 'vy', 'valpha');

if is_octave()
    waitforbuttonpress;
end

% check consistency
xp = x_sim(1,:);
yp = x_sim(2,:);
check = xp.^2 + yp.^2 - length_pendulum^2;
if any( max(abs(check)) > 1e-10 )
    disp('note: check for constant pendulum length failed');
end

