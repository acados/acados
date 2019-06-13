%% test of native matlab interface
clear all



% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	disp('ERROR: env.sh has not been sourced! Before executing this example, run:');
	disp('source env.sh');
	return;
end



%% arguments
compile_mex = 'true';
codgen_model = 'true';
method = 'erk';
%method = 'irk';
sens_forw = 'false';
num_stages = 4;
num_steps = 4;

h = 0.1;
x0 = [0; 1e-1; 0; 0e0];
u = 0;



%% model
model = pendulum_on_cart_model;

nx = model.nx;
nu = model.nu;



%% acados sim model
sim_model = acados_sim_model();
sim_model.set('T', h);
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

%sim_model.model_struct



%% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('compile_mex', compile_mex);
sim_opts.set('codgen_model', codgen_model);
sim_opts.set('num_stages', num_stages);
sim_opts.set('num_steps', num_steps);
sim_opts.set('method', method);
sim_opts.set('sens_forw', sens_forw);

%sim_opts.opts_struct



%% acados sim
% create sim
sim = acados_sim(sim_model, sim_opts);
% (re)set numerical part of model
%sim.set('T', 0.5);
%sim.C_sim
%sim.C_sim_ext_fun


N_sim = 100;

x_sim = zeros(nx, N_sim+1);
x_sim(:,1) = x0;

tic
for ii=1:N_sim
	
	% set initial state
	sim.set('x', x_sim(:,ii));
	sim.set('u', u);

	% solve
	sim.solve();


	% get simulated state
	x_sim(:,ii+1) = sim.get('xn');

end
simulation_time = toc


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

%x_sim

for ii=1:N_sim+1
	x_cur = x_sim(:,ii);
	visualize;
end

figure(2);
plot(1:N_sim+1, x_sim);
legend('p', 'theta', 'v', 'omega');


fprintf('\nsuccess!\n\n');


waitforbuttonpress;


return;

