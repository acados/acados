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



% NOTE: `acados` currently supports both an old MATLAB/Octave interface (< v0.4.0)
% as well as a new interface (>= v0.4.0).

% THIS EXAMPLE still uses the OLD interface. If you are new to `acados` please start
% with the examples that have been ported to the new interface already.
% see https://github.com/acados/acados/issues/1196#issuecomment-2311822122)


clear all; clc;
check_acados_requirements();

%% arguments
compile_interface = 'auto';
%method = 'erk';
method = 'irk';
method = 'irk_gnsf';
gnsf_detect_struct = 'true';
sens_forw = 'true';
num_stages = 4;
num_steps = 4;

%% model
model = linear_mass_spring_model();

nx = model.nx;
nu = model.nu;

%% acados sim model
sim_model = acados_sim_model();
sim_model.set('T', 0.5);
if (strcmp(method, 'erk'))
	sim_model.set('dyn_type', 'explicit');
	sim_model.set('dyn_expr_f', model.dyn_expr_f_expl);
	sim_model.set('sym_x', model.sym_x);
	if isfield(model, 'sym_u')
		sim_model.set('sym_u', model.sym_u);
	end
else % irk irk_gnsf
	sim_model.set('dyn_type', 'implicit');
	sim_model.set('dyn_expr_f', model.dyn_expr_f_impl);
	sim_model.set('sym_x', model.sym_x);
	sim_model.set('sym_xdot', model.sym_xdot);
	if isfield(model, 'sym_u')
		sim_model.set('sym_u', model.sym_u);
	end
end

%% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('compile_interface', compile_interface);
sim_opts.set('num_stages', num_stages);
sim_opts.set('num_steps', num_steps);
sim_opts.set('method', method);
sim_opts.set('sens_forw', sens_forw);
if (strcmp(method, 'irk_gnsf'))
	sim_opts.set('gnsf_detect_struct', gnsf_detect_struct);
end

%% acados sim
% create sim
sim_solver = acados_sim(sim_model, sim_opts);
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

