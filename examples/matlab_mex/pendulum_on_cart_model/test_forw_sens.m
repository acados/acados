%% test of native matlab interface
clear all



%% arguments
compile_mex = 'true';
codgen_model = 'true';
method = 'irk';
sens_forw = 'true';
num_stages = 4;
num_steps = 4;

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

% solve
sim.solve();

% xn
xn = sim.get('xn');
xn
% S_forw
S_forw_ind = sim.get('S_forw');
S_forw_ind



%% compute forward sensitivities using finite differences

% TODO set sens_forw to 0

S_forw_fd = zeros(nx, nx+nu);

% asymmetric finite differences

for ii=1:nx

	dx = zeros(nx, 1);
	dx(ii) = 1.0;

	sim.set('x', x0+epsilon*dx);
	sim.set('u', u);

	% solve
	sim.solve();

	% xn
	xn_tmp = sim.get('xn');

	S_forw_fd(:,ii) = (xn_tmp - xn) / epsilon;

end

for ii=1:nu

	du = zeros(nu, 1);
	du(ii) = 1.0;

	sim.set('x', x0);
	sim.set('u', u+epsilon*du);

	% solve
	sim.solve();

	% xn
	xn_tmp = sim.get('xn');

	S_forw_fd(:,nx+ii) = (xn_tmp - xn) / epsilon;

end

S_forw_fd


%% compute error

error_abs = S_forw_fd - S_forw_ind



fprintf('\nsuccess!\n\n');


return;


