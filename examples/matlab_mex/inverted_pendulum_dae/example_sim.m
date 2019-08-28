%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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

%% test of native matlab interface
clear VARIABLES

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	disp('ERROR: env.sh has not been sourced! Before executing this example, run:');
	disp('source env.sh');
	return;
end

%% arguments
compile_mex = 'true'; % true, false
codgen_model = 'true'; % true, false
gnsf_detect_struct = 'true'; % true, false
method = 'irk'; % irk, irk_gnsf, [erk]
sens_forw = 'true'; % true, false
jac_reuse = 'false'; % true, false
num_stages = 3;
num_steps = 3;
newton_iter = 3;
model_name = 'inv_pend_dae';

x0 = [1; -5; 1; 0.1; -0.5; 0.1];
u = 1;

%% model
model = inverted_pendulum_dae_model;
disp('state')
disp(model.sym_x)

nx = length(model.sym_x);
nu = length(model.sym_u);
nz = length(model.sym_z);

%% acados sim model
sim_model = acados_sim_model();
sim_model.set('name', model_name);
sim_model.set('T', 0.01); % simulation time

sim_model.set('sym_x', model.sym_x);
if isfield(model, 'sym_u')
    sim_model.set('sym_u', model.sym_u);
end
if isfield(model, 'sym_p')
    sim_model.set('sym_p', model.sym_p);
end
sim_model.set('dim_nx', nx);
sim_model.set('dim_nu', nu);


% Note: DAEs can only be used with implicit integrator
sim_model.set('dyn_type', 'implicit');
sim_model.set('dyn_expr_f', model.expr_f_impl);
sim_model.set('sym_xdot', model.sym_xdot);
if isfield(model, 'sym_z')
	sim_model.set('sym_z', model.sym_z);
end
sim_model.set('dim_nz', nz);

%% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('compile_mex', compile_mex);
sim_opts.set('codgen_model', codgen_model);
sim_opts.set('num_stages', num_stages);
sim_opts.set('num_steps', num_steps);
sim_opts.set('newton_iter', newton_iter);
sim_opts.set('method', method);
sim_opts.set('sens_forw', sens_forw);
sim_opts.set('sens_adj', 'true');
sim_opts.set('sens_algebraic', 'true');
sim_opts.set('output_z', 'true');
sim_opts.set('sens_hess', 'false');
sim_opts.set('jac_reuse', jac_reuse);
if (strcmp(method, 'irk_gnsf'))
	sim_opts.set('gnsf_detect_struct', gnsf_detect_struct);
end


%% acados sim
% create integrator
sim = acados_sim(sim_model, sim_opts);

N_sim = 1000;

x_sim = zeros(nx, N_sim+1);
x_sim(:,1) = x0;

tic
for ii=1:N_sim
	
	% set initial state
	sim.set('x', x_sim(:,ii));
	sim.set('u', u);

	% set adjoint seed
	sim.set('seed_adj', ones(nx,1));

    % initialize implicit integrator
    if (strcmp(method, 'irk'))
        sim.set('xdot', zeros(nx,1));
        sim.set('z', zeros(nz,1));
    elseif (strcmp(method, 'irk_gnsf'))
        n_out = sim.model_struct.dim_gnsf_nout;
        sim.set('phi_guess', zeros(n_out,1));
    end

	% solve
	sim.solve();

	% get simulated state
	x_sim(:,ii+1) = sim.get('xn');

	% forward sensitivities
	S_forw = sim.get('S_forw');
	Sx = sim.get('Sx');
	Su = sim.get('Su');
	
end
format short e
xfinal = sim.get('xn')'
S_forw
S_adj = sim.get('S_adj')'
z = sim.get('zn')' % approximate value of algebraic variables at start of simulation
S_alg = sim.get('S_algebraic') % sensitivities of algebraic variables z

simulation_time = toc


figure(2);
plot(1:N_sim+1, x_sim);
legend('xpos', 'ypos', 'alpha', 'vx', 'vy', 'valpha');

waitforbuttonpress;

