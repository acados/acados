%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

%% test of native matlab interface
clear VARIABLES

addpath('../linear_mass_spring_model/');

%% arguments
compile_interface = 'auto';
method = 'irk';
sens_forw = 'true';
num_stages = 4;
num_steps = 4;

Ts = 0.1;

%% model
model = linear_mass_spring_model();

model_name = ['lin_mass_' method];
nx = model.nx;
nu = model.nu;

%% acados sim model
sim_model = acados_sim_model();
sim_model.set('T', Ts);
sim_model.set('name', model_name);

sim_model.set('sym_x', model.sym_x);
if isfield(model, 'sym_u')
    sim_model.set('sym_u', model.sym_u);
end
if isfield(model, 'sym_p')
    sim_model.set('sym_p', model.sym_p);
end


if (strcmp(method, 'erk'))
    sim_model.set('dyn_type', 'explicit');
    sim_model.set('dyn_expr_f', model.dyn_expr_f_expl);
else % irk irk_gnsf
    sim_model.set('dyn_type', 'implicit');
    sim_model.set('dyn_expr_f', model.dyn_expr_f_impl);
    sim_model.set('sym_xdot', model.sym_xdot);
end


%% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('compile_interface', compile_interface);
sim_opts.set('num_stages', num_stages);
sim_opts.set('num_steps', num_steps);
sim_opts.set('method', method);
sim_opts.set('sens_forw', sens_forw);

%% acados sim
% create sim
sim_solver = acados_sim(sim_model, sim_opts);

% Note: this does not work with gnsf, because it needs to be available
% in the precomputation phase
% 	sim_solver.set('T', Ts);

%% test check, this should fail!
% set initial state
sim_solver.set('x', zeros(nx+1, 1));
