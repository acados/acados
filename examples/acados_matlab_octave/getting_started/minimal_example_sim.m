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

%% minimal example of acados integrator matlab interface
clear all

addpath('../pendulum_on_cart_model')

check_acados_requirements()

%% arguments
compile_interface = 'auto';
method = 'irk_gnsf'; % irk, irk_gnsf
model_name = 'sim_pendulum';

% simulation parameters
N_sim = 100;
h = 0.1; % simulation time
x0 = [0; 1e-1; 0; 0]; % initial state
u0 = 0; % control input

%% define model dynamics
model = pendulum_on_cart_model();

nx = model.nx;
nu = model.nu;

%% acados sim model
sim_model = acados_sim_model();
sim_model.set('name', model_name);
sim_model.set('T', h);

sim_model.set('sym_x', model.sym_x);
if isfield(model, 'sym_u')
    sim_model.set('sym_u', model.sym_u);
end

% explit integrator (erk) take explicit ODE expression
if (strcmp(method, 'erk'))
	sim_model.set('dyn_type', 'explicit');
	sim_model.set('dyn_expr_f', model.expr_f_expl);
else % implicit integrators (irk irk_gnsf) take implicit ODE expression
	sim_model.set('dyn_type', 'implicit');
	sim_model.set('dyn_expr_f', model.expr_f_impl);
	sim_model.set('sym_xdot', model.sym_xdot);
end

%% acados sim opts
sim_opts = acados_sim_opts();

sim_opts.set('compile_interface', compile_interface);
sim_opts.set('num_stages', 2);
sim_opts.set('num_steps', 3);
sim_opts.set('newton_iter', 2); % for implicit intgrators
sim_opts.set('method', method);
sim_opts.set('sens_forw', 'true'); % generate forward sensitivities
if (strcmp(method, 'irk_gnsf'))
	sim_opts.set('gnsf_detect_struct', 'true');
end


%% create integrator
sim = acados_sim(sim_model, sim_opts);


%% simulate system in loop
x_sim = zeros(nx, N_sim+1);
x_sim(:,1) = x0;

for ii=1:N_sim
	
	% set initial state
	sim.set('x', x_sim(:,ii));
	sim.set('u', u0);

    % initialize implicit integrator
    if (strcmp(method, 'irk'))
        sim.set('xdot', zeros(nx,1));
    elseif (strcmp(method, 'irk_gnsf'))
        n_out = sim.model_struct.dim_gnsf_nout;
        sim.set('phi_guess', zeros(n_out,1));
    end

	% solve
	sim.solve();

	% get simulated state
	x_sim(:,ii+1) = sim.get('xn');
end

% forward sensitivities ( dxn_d[x0,u] )
S_forw = sim.get('S_forw');

%% plots
% for ii=1:N_sim+1
% 	x_cur = x_sim(:,ii);
% 	visualize;
% end

figure;
plot(1:N_sim+1, x_sim);
legend('p', 'theta', 'v', 'omega');

%% dump gnsf
% sim.model_struct.name = 'pendulum_ode';
% dump_gnsf_functions(sim.model_struct)
