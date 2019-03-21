%% test of native matlab interface
clear all



%% arguments
compile_mex = 'true';
codgen_model = 'true';
method = 'irk';
sens_forw = 'true';
num_stages = 4;
num_steps = 4;



%% model
model = linear_mass_spring_model;

nx = model.nx;
nu = model.nu;



%% acados sim model
sim_model = acados_sim_model();
sim_model.set('T', 0.5);
if (strcmp(method, 'erk'))
	sim_model.set('dyn_type', 'explicit');
	sim_model.set('dyn_expr_f', model.expr_f_expl);
	sim_model.set('sym_x', model.sym_x);
	if isfield(model, 'sym_u')
		sim_model.set('sym_u', model.sym_u);
	end
	sim_model.set('dim_nx', model.nx);
	sim_model.set('dim_nu', model.nu);
else % irk
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
	sim_model.set('dim_nx', model.nx);
	sim_model.set('dim_nu', model.nu);
%	sim_model.set('nz', model.nz);
end

sim_model.model_struct




%% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('compile_mex', compile_mex);
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
sim.C_sim
sim.C_sim_ext_fun



x0 = ones(nx, 1); %x0(1) = 2.0;
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
