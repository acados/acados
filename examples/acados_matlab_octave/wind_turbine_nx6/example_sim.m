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

%% test of native matlab interface
clear all


% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end

% load sim data
load testSim.mat


%% arguments
compile_interface = 'auto';
codgen_model = 'true';
method = 'irk'; % irk, erk, irk_gnsf
% method = 'irk_gnsf';
sens_forw = 'true';
num_stages = 4;
num_steps = 4;

%% parametric model
model = sim_model_wind_turbine_nx6;
nx = model.nx;
nu = model.nu;
np = model.np;

%% acados sim model
Ts = 0.2;
sim_model = acados_sim_model();
sim_model.set('T', Ts);

sim_model.set('sym_x', model.sym_x);
if isfield(model, 'sym_u')
    sim_model.set('sym_u', model.sym_u);
end
if isfield(model, 'sym_p')
    sim_model.set('sym_p', model.sym_p);
end

if (strcmp(method, 'erk'))
	sim_model.set('dyn_type', 'explicit');
	sim_model.set('dyn_expr_f', model.expr_f_expl);
else % irk irk_gnsf
	sim_model.set('dyn_type', 'implicit');
	sim_model.set('dyn_expr_f', model.expr_f_impl);
	sim_model.set('sym_xdot', model.sym_xdot);

%	if isfield(model, 'sym_z')
%		sim_model.set('sym_z', model.sym_z);
%	end
end

%% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('compile_interface', compile_interface);
sim_opts.set('codgen_model', codgen_model);
sim_opts.set('num_stages', num_stages);
sim_opts.set('num_steps', num_steps);
sim_opts.set('method', method);
sim_opts.set('sens_forw', sens_forw);


%% acados sim
% create sim
sim = acados_sim(sim_model, sim_opts);

% to avoid unstable behavior introduce a small pi-contorller for rotor speed tracking
uctrl = 0.0;
uctrlI = 0.0;
kI = 1e-1;
kP = 10;


nsim = 15;

x_sim = zeros(nx, nsim+1);
x_sim(:,1) = statesFAST(1,:);

tic;
for nn=1:nsim

	% compute input
	u = Usim(nn,1:2);
	u(2) = max(u(2) - uctrl, 0);

	% update state, input, parameter
	sim.set('x', x_sim(:,nn));
	sim.set('u', u);
	sim.set('p', Usim(nn,3));

	% solve
	sim.solve();

	x_sim(:,nn+1) = sim.get('xn');

	% update PI contoller
	ctrlErr = statesFAST(nn+1,1) - x_sim(1,nn+1);
	uctrlI = uctrlI + kI*ctrlErr*Ts;
	uctrl = kP*ctrlErr + uctrlI;

end

time_solve = toc/nsim

%statesFAST(1:nsim+1,:)'
x_sim(:,1:nsim+1)

%S_forw = sim.get('S_forw');
%Sx = sim.get('Sx');
%Su = sim.get('Su');


fprintf('\nsuccess!\n\n');

