%% test of native matlab interface
clear all



% load sim data
load testSim.mat



%% arguments
compile_mex = 'true';
codgen_model = 'true';
method = 'irk';
sens_forw = 'true';
num_stages = 4;
num_steps = 4;



%% parametric model
model = sim_model_wind_turbine_nx6;

model

nx = model.nx;
nu = model.nu;
np = model.np;



%% acados sim model
Ts = 0.2;
sim_model = acados_sim_model();
sim_model.set('T', Ts);
if (strcmp(method, 'erk'))
	sim_model.set('dyn_type', 'explicit');
	sim_model.set('dyn_param_f', 'true');
	sim_model.set('dyn_expr_f', model.expr_f_expl);
	sim_model.set('sym_x', model.sym_x);
	if isfield(model, 'sym_u')
		sim_model.set('sym_u', model.sym_u);
	end
	if isfield(model, 'sym_p')
		sim_model.set('sym_p', model.sym_p);
	end
	sim_model.set('dim_nx', model.nx);
	sim_model.set('dim_nu', model.nu);
	sim_model.set('dim_np', model.np);
else % irk
	sim_model.set('dyn_type', 'implicit');
	sim_model.set('dyn_param_f', 'true');
	sim_model.set('dyn_expr_f', model.expr_f_impl);
	sim_model.set('sym_x', model.sym_x);
	sim_model.set('sym_xdot', model.sym_xdot);
	if isfield(model, 'sym_u')
		sim_model.set('sym_u', model.sym_u);
	end
%	if isfield(model, 'sym_z')
%		sim_model.set('sym_z', model.sym_z);
%	end
	if isfield(model, 'sym_p')
		sim_model.set('sym_p', model.sym_p);
	end
	sim_model.set('dim_nx', model.nx);
	sim_model.set('dim_nu', model.nu);
%	sim_model.set('dim_nz', model.nz);
	sim_model.set('dim_np', model.np);
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

% get TODO with return value !!!!!
% xn
%xn = sim.get('xn');
%xn
% S_forw
%S_forw = sim.get('S_forw');
%S_forw
% Sx
%Sx = sim.get('Sx');
%Sx
% Su
%Su = sim.get('Su');
%Su


fprintf('\nsuccess!\n\n');


return;
