%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

addpath('../pendulum_on_cart_model/');

for integrator = {'ERK', 'IRK'}
	method = integrator{1};
	num_stages = 4;
	num_steps = 3;

	Ts = 0.1;
	x0 = [1e-1; 1e0; 2e-1; 2e0];
	u = 0;
	FD_epsilon = 1e-6;

	model = get_pendulum_on_cart_model();

	model_name = ['pendulum_' method];
	nx = length(model.x);
	nu = length(model.u);

	sim = AcadosSim();
	sim.model = model;
	sim.model.name = model_name;
	sim.solver_options.Tsim = Ts;
	sim.solver_options.integrator_type = method;
	sim.solver_options.num_stages = num_stages;
	sim.solver_options.num_steps = num_steps;
	sim.solver_options.sens_forw = true;
	sim.solver_options.sens_adj = true;
	sim.solver_options.sens_hess = true;
	if strcmp(method, 'ERK')
	    sim.model.f_expl_expr = model.f_expl_expr;
	else
	    sim.model.f_impl_expr = model.f_impl_expr;
	end
	sim_solver = AcadosSimSolver(sim);

	S_hess_ind = zeros(nx+nu, nx+nu, nx);

	S_hess_fd = zeros(nx+nu, nx+nu, nx);

	for jj=1:nx
		sim_solver.set('x', x0);
		sim_solver.set('u', u);

		lambda = zeros(nx, 1);
		lambda(jj) = 1.0;
		sim_solver.set('seed_adj', lambda);

		sim_solver.solve();

		S_hess = sim_solver.get('S_hess');
		S_hess_ind(:, :, jj) = S_hess;

		S_adj = sim_solver.get('S_adj');

		for ii=1:nx
			dx = zeros(nx, 1);
			dx(ii) = 1.0;

			sim_solver.set('x', x0+FD_epsilon*dx);
			sim_solver.set('u', u);

			sim_solver.solve();
			S_adj_tmp = sim_solver.get('S_adj');
			S_hess_fd(:, ii, jj) = (S_adj_tmp - S_adj) / FD_epsilon;
		
		end

		for ii=1:nu
			du = zeros(nu, 1);
			du(ii) = 1.0;

			sim_solver.set('x', x0);
			sim_solver.set('u', u+FD_epsilon*du);

			sim_solver.solve();
			S_adj_tmp = sim_solver.get('S_adj');
			S_hess_fd(:, nx+ii, jj) = (S_adj_tmp - S_adj) / FD_epsilon;
        end
	end

	error_abs = max(max(max(abs(S_hess_fd - S_hess_ind))));
	disp(' ')
	disp(['integrator:  ' method]);
	disp(['error hessian (wrt finite differences):   ' num2str(error_abs)])
	disp(' ')
    if error_abs > 1e-6
        error(strcat('test_sens_hess FAIL: second order sensitivity error too large: \n',...
			'for integrator:\t', method));
	end
end

fprintf('\nTEST_HESSIANS: success!\n\n');

return;
