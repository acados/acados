%% test of native matlab interface
clear all



%% arguments
compile_mex = 'true';
codgen_model = 'true';
method = 'erk';
sens_forw = 'true';
sens_adj = 'true';
sens_hess = 'true';
num_stages = 4;
num_steps = 3;

Ts = 0.1;
x0 = [1e-1; 1e0; 2e-1; 2e0];
u = 0;
epsilon = 1e-6;



%% model
model = pendulum_on_cart_model;

nx = model.nx;
nu = model.nu;



%% acados sim model
sim_model = acados_sim_model();
%sim_model.set('T', Ts);
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
sim_opts.set('sens_adj', sens_adj);
sim_opts.set('sens_hess', sens_hess);

%sim_opts.opts_struct



%% acados sim
% create sim
sim = acados_sim(sim_model, sim_opts);
% (re)set numerical part of model
%sim.C_sim
%sim.C_sim_ext_fun



% set simulation time
sim.set('T', Ts);

% set initial state
sim.set('x', x0);
sim.set('u', u);



%% compute hessian sensitivities using internal numerical differentiation

S_hess_ind = zeros(nx+nu, nx+nu, nx);

% asymmetric finite differences

for ii=1:nx

	dx = zeros(nx, 1);
	dx(ii) = 1.0;

	sim.set('seed_adj', dx);

	% solve
	sim.solve();

	% xn
	S_hess = sim.get('S_hess');

	S_hess_ind(:, :, ii) = S_hess;

end

S_hess_ind


%% compute error

%error_abs = S_forw_adj - S_forw_ind



fprintf('\nsuccess!\n\n');


return;




