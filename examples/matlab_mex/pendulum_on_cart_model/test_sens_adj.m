%% test of native matlab interface
clear VARIABLES

for integrator = {'irk_gnsf', 'irk', 'erk'}
	%% arguments
	compile_mex = 'true';
	codgen_model = 'true';
	method = integrator{1}; %'irk'; 'irk_gnsf'; 'erk';
	sens_forw = 'true';
	sens_adj = 'true';
	num_stages = 4;
	num_steps = 3;

	Ts = 0.1;
	x0 = [1e-1; 1e0; 2e-1; 2e0];
	u = 0;

    disp(['testing integrator  ' method])
	%% model
	model = pendulum_on_cart_model;

	model_name = ['pendulum_' method];
	nx = model.nx;
	nu = model.nu;

	%% acados sim model
	sim_model = acados_sim_model();
	sim_model.set('T', Ts);
	sim_model.set('name', model_name);
	if (strcmp(method, 'erk'))
		sim_model.set('dyn_type', 'explicit');
		sim_model.set('dyn_expr_f', model.expr_f_expl);
		sim_model.set('sym_x', model.sym_x);
		if isfield(model, 'sym_u')
			sim_model.set('sym_u', model.sym_u);
		end
		sim_model.set('dim_nx', model.nx);
		sim_model.set('dim_nu', model.nu);
	elseif (strcmp(method, 'irk')||strcmp(method, 'irk_gnsf'))
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

	%% acados sim opts
	sim_opts = acados_sim_opts();
	sim_opts.set('compile_mex', compile_mex);
	sim_opts.set('codgen_model', codgen_model);
	sim_opts.set('num_stages', num_stages);
	sim_opts.set('num_steps', num_steps);
	sim_opts.set('method', method);
	sim_opts.set('sens_forw', sens_forw);
	sim_opts.set('sens_adj', sens_adj);


	%% acados sim
	% create sim
	sim = acados_sim(sim_model, sim_opts);

    % Note: this does not work with gnsf, because it needs to be available
    % in the precomputation phase
    % 	sim.set('T', Ts);

	% set initial state
	sim.set('x', x0);
	sim.set('u', u);

	% solve
	sim.solve();

	xn = sim.get('xn');
	S_forw_ind = sim.get('S_forw');

	%% compute forward sensitivities through adjoint sensitivities with unit seeds
	sim_opts.set('sens_forw', 'false');
	sim_opts.set('sens_adj', 'true');

	S_forw_adj = zeros(nx, nx+nu);
	for ii=1:nx
		% set seed
		dx = zeros(nx, 1);
		dx(ii) = 1.0;
		sim.set('seed_adj', dx);

		sim.solve();

		S_adj = sim.get('S_adj');
		S_forw_adj(ii,:) = S_adj;
	end

	%% compute & check error
	error_abs = max(max(max(abs(S_forw_adj - S_forw_ind))));
	disp(' ')
	disp(['integrator:  ' method]);
	disp(['error adjoint sensitivities (wrt forward sens):   ' num2str(error_abs)])
	disp(' ')
	if error_abs > 1e-14
		disp(['adjoint sensitivities error too large -> debug mode'])
		keyboard
	end
end

fprintf('\nTEST_ADJOINT_SENSITIVITIES: success!\n\n');

return;