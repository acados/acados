%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

clear VARIABLES

addpath('../linear_mass_spring_model/');

for integrator = {'GNSF', 'IRK', 'ERK'}
    method = integrator{1};

    sens_forw = true;
    jac_reuse = true;
    num_stages = 3;
    num_steps = 4;
    newton_iter = 3;

    Ts = 0.1;
    FD_epsilon = 1e-6;

    model = get_linear_mass_spring_model();

    model_name = ['lin_mass_' method];
    nx = length(model.x);
    nu = length(model.u);
    x0 = ones(nx,1);
    u = ones(nu,1);

    sim = AcadosSim();
    sim.model = model;
    sim.model.name = model_name;
    sim.solver_options.Tsim = Ts;
    sim.solver_options.integrator_type = method;
    sim.solver_options.num_stages = num_stages;
    sim.solver_options.num_steps = num_steps;
    sim.solver_options.newton_iter = newton_iter;
    sim.solver_options.sens_forw = sens_forw;
    sim.solver_options.jac_reuse = jac_reuse;
    if strcmp(method, 'ERK')
        sim.model.f_expl_expr = model.f_expl_expr;
    else
        sim.model.f_impl_expr = model.f_impl_expr;
    end
    sim_solver = AcadosSimSolver(sim);

    sim_solver.set('x', x0);
    sim_solver.set('u', u);

    if strcmp(method, 'IRK')
        sim_solver.set('xdot', zeros(nx,1));
    elseif strcmp(method, 'GNSF')
        n_out = sim_solver.sim.model.gnsf_model.dims.nout;
        sim_solver.set('phi_guess', zeros(n_out,1));
    end

    sim_solver.solve();

    xn = sim_solver.get('xn');
    S_forw_ind = sim_solver.get('S_forw');

    S_forw_fd = zeros(nx, nx+nu);

    for ii=1:nx
        dx = zeros(nx, 1);
        dx(ii) = 1.0;

        sim_solver.set('x', x0+FD_epsilon*dx);
        sim_solver.set('u', u);
        sim_solver.solve();

        xn_tmp = sim_solver.get('xn');
        S_forw_fd(:,ii) = (xn_tmp - xn) / FD_epsilon;
    end

    for ii=1:nu
        du = zeros(nu, 1);
        du(ii) = 1.0;

        sim_solver.set('x', x0);
        sim_solver.set('u', u+FD_epsilon*du);
        sim_solver.solve();

        xn_tmp = sim_solver.get('xn');
        S_forw_fd(:,nx+ii) = (xn_tmp - xn) / FD_epsilon;
    end

    error_abs = max(max(abs(S_forw_fd - S_forw_ind)));
    disp(' ')
    disp(['integrator:  ' method]);
    disp(['error forward sensitivities (wrt finite differences):   ' num2str(error_abs)])
    disp(' ')
    if error_abs > 1e-6
        error(strcat('test_sens_forw FAIL: forward sensitivities error too large: \n',...
            'for integrator:\t', method));
    end
end

fprintf('\nTEST_FORWARD_SENSITIVITIES: success!\n\n');
