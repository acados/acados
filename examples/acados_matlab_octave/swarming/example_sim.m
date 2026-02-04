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

% Author: Enrica


% NOTE: `acados` currently supports both an old MATLAB/Octave interface (< v0.4.0)
% as well as a new interface (>= v0.4.0).

% THIS EXAMPLE still uses the OLD interface. If you are new to `acados` please start
% with the examples that have been ported to the new interface already.
% see https://github.com/acados/acados/issues/1196#issuecomment-2311822122)



% This file allows the simulation of the dynamics of a swarm of robots.
% Here, the swarm is composed by N agents with decoupled, linear dynamics.
%

%% Test of native matlab interface
clear all

% Check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	disp('ERROR: env.sh has not been sourced! Before executing this example, run:');
	disp('source env.sh');
	return;
end

%% Arguments

% Time parameters
dt = 0.1; % discretization step
T = 10; % total time of the simulation
nb_steps = floor(T/dt); % number of time steps during the simulation

% Structure S with the swarming parameters
S.N = 3; % number of agents in the swarm
S.d_ref = 5; % reference distance among every couple of neighboring agents
S.u_ref = [1;0;0]; % reference direction of velocity for all agents
S.v_ref = 6; % reference speed for all agents

% Rename parameters
N = S.N;

% Initial conditions
x0 = [10*rand(3*N,1); 2*rand(3*N,1)]; % initial condition (3D positions ...
    % and velocities of N agents)
u = zeros(N*3,1); % control input is null (acceleration)

% Simulation parameters
compile_interface = 'auto';
%gnsf_detect_struct = 'true';
method = 'erk';
%method = 'irk';
%method = 'irk_gnsf';
sens_forw = 'true';
num_stages = 4;
num_steps = 4;
model_name = 'sim_swarming';

%% Model

model = swarming_model(S);

nx = model.nx;
nu = model.nu;

%% Acados simutation model

sim_model = acados_sim_model();
sim_model.set('name', model_name);
sim_model.set('T', dt);
if (strcmp(method, 'erk'))
	sim_model.set('dyn_type', 'explicit');
	sim_model.set('dyn_expr_f', model.expr_f_expl);
	sim_model.set('sym_x', model.sym_x);
	if isfield(model, 'sym_u')
		sim_model.set('sym_u', model.sym_u);
	end
else % irk irk_gnsf
	sim_model.set('dyn_type', 'implicit');
	sim_model.set('dyn_expr_f', model.expr_f_impl);
	sim_model.set('sym_x', model.sym_x);
	sim_model.set('sym_xdot', model.sym_xdot);
	if isfield(model, 'sym_u')
		sim_model.set('sym_u', model.sym_u);
	end
%	if isfield(model, 'sym_z')
%		sim_model.set('sym_z', model.sym_z);
%	end
end

%% Acados simutation options
sim_opts = acados_sim_opts();
sim_opts.set('compile_interface', compile_interface);
sim_opts.set('num_stages', num_stages);
sim_opts.set('num_steps', num_steps);
sim_opts.set('method', method);
sim_opts.set('sens_forw', sens_forw);
if (strcmp(method, 'irk_gnsf'))
	sim_opts.set('gnsf_detect_struct', gnsf_detect_struct);
end

%% Acados simulation

% Create sim
sim_solver = acados_sim(sim_model, sim_opts);
% (Re)set numerical part of model
%sim_solver.set('T', 0.5);
%sim_solver.C_sim
%sim_solver.C_sim_ext_fun

x_history = zeros(nb_steps+1, nx);
x_history(1,:) = x0';

tic
for k = 1:nb_steps
	% Set initial state
	sim_solver.set('x', x_history(k,:));
	sim_solver.set('u', u);

	% Solve
	sim_solver.solve();

	% Get simulated state
	x_history(k+1,:) = sim_solver.get('xn');
end
simulation_time = toc;
disp(strcat('Simulation time: ',num2str(simulation_time)));

xn = sim_solver.get('xn');
S_forw = sim_solver.get('S_forw');
Sx = sim_solver.get('Sx');
Su = sim_solver.get('Su');

%% Plot variables

fontsize = 12; % font size for plots

% Plot
pos_history = x_history(:,1:3*N);
vel_history = x_history(:,(3*N+1):end);
for agent = 1:N
    plot3(pos_history(:,(agent-1)*3+2), pos_history(:,(agent-1)*3+1), ...
        - pos_history(:,(agent-1)*3+3));
    hold on;
end
% title('Swarm trajectories');
xlabel('Y Position [m]','fontsize',fontsize);
ylabel('X Position [m]','fontsize',fontsize);
zlabel('Z Position [m]','fontsize',fontsize);
view(2);

%% End of simulation

fprintf('\nsuccess!\n\n');
return;

