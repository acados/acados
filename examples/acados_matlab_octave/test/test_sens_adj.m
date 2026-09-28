%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

clear VARIABLES

addpath('../pendulum_on_cart_model/');

for integrator = {'GNSF', 'IRK', 'ERK'}
    method = integrator{1};
    num_stages = 4;
    num_steps = 4;

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
	if strcmp(method, 'ERK')
	    sim.model.f_expl_expr = model.f_expl_expr;
	else
	    sim.model.f_impl_expr = model.f_impl_expr;
	end
	sim_solver = AcadosSimSolver(sim);

	sim_solver.set('x', x0);
	sim_solver.set('u', u);

	sim_solver.solve();

	xn = sim_solver.get('xn');
	S_forw_ind = sim_solver.get('S_forw');

	S_forw_adj = zeros(nx, nx+nu);
	for ii=1:nx
		dx = zeros(nx, 1);
		dx(ii) = 1.0;
		sim_solver.set('seed_adj', dx);

		sim_solver.solve();

		S_adj = sim_solver.get('S_adj');
		S_forw_adj(ii,:) = S_adj;
	end

	error_abs = max(max(max(abs(S_forw_adj - S_forw_ind))));
	disp(' ')
	disp(['integrator:  ' method]);
	disp(['error adjoint sensitivities (wrt forward sens):   ' num2str(error_abs)])
	disp(' ')
	if error_abs > 1e-14
        error(strcat('test_sens_adj FAIL: adjoint sensitivities error too large: \n',...
            'for integrator:\t', method));
	end
end

fprintf('\nTEST_ADJOINT_SENSITIVITIES: success!\n\n');
