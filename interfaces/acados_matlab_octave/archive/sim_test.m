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



%% arguments
compile_interface = 'true';
codgen_model = 'true';
method = 'irk';
sens_forw = 'true';
num_stages = 4;
num_steps = 4;



%% model
%model = linear_model;
model = linear_mass_spring_model;
%model = crane_model;

nx = model.nx;
nu = model.nu;



%% acados sim model
sim_model = acados_sim_model();
sim_model.set('T', 0.5);
if (strcmp(method, 'erk'))
	sim_model.set('dyn_type', 'explicit');
	sim_model.set('expr_f', model.expr_f_expl);
	sim_model.set('sym_x', model.sym_x);
	if isfield(model, 'sym_u')
		sim_model.set('sym_u', model.sym_u);
	end
	sim_model.set('nx', model.nx);
	sim_model.set('nu', model.nu);
else % irk
	sim_model.set('dyn_type', 'implicit');
	sim_model.set('expr_f', model.expr_f_impl);
	sim_model.set('sym_x', model.sym_x);
	sim_model.set('sym_xdot', model.sym_xdot);
	if isfield(model, 'sym_u')
		sim_model.set('sym_u', model.sym_u);
	end
%	if isfield(model, 'sym_z')
%		sim_model.set('sym_z', model.sym_z);
%	end
	sim_model.set('nx', model.nx);
	sim_model.set('nu', model.nu);
%	sim_model.set('dim_nz', model.nz);
end

sim_model.model_struct




%% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('compile_interface', compile_interface);
sim_opts.set('codgen_model', codgen_model);
sim_opts.set('num_stages', num_stages);
sim_opts.set('num_steps', num_steps);
sim_opts.set('method', method);
sim_opts.set('sens_forw', sens_forw);

sim_opts.opts_struct



%% acados sim
% create sim
sim = acados_sim(sim_model, sim_opts);
% (re)set numerical part of model
%sim.set('T', 0.5);

x0 = ones(sim_model.nx, 1); %x0(1) = 2.0;
tic;
sim.set('x', x0);
time_set_x = toc

u = ones(nu, 1);
sim.set('u', u);

% solve
tic;
sim.solve();
time_solve = toc


% get TODO with return value !!!!!
% xn
xn = sim.get('xn');
xn
% S_forw
S_forw = sim.get('S_forw');
S_forw
% Sx
Sx = sim.get('Sx');
Sx
% Su
Su = sim.get('Su');
Su


fprintf('\nsuccess!\n\n');


return;
