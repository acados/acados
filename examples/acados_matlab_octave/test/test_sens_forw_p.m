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

    sens_forw   = true;
    jac_reuse   = true;
    num_stages  = 3;
    num_steps   = 4;
    newton_iter = 3;

    Ts         = 0.1;
    x0         = [1e-1; 1e0; 2e-1; 2e0];
    u          = 0;
    FD_epsilon = 1e-6;

    old_model = pendulum_on_cart_model_with_param();
    model = AcadosModel();
    model.name = ['pendulum_sens_p_' method];
    model.x = old_model.sym_x;
    model.xdot = old_model.sym_xdot;
    model.u = old_model.sym_u;
    model.p = old_model.sym_p;
    model.f_expl_expr = old_model.dyn_expr_f_expl;
    model.f_impl_expr = old_model.dyn_expr_f_impl;
    model_name = ['pendulum_sens_p_' method];

    nx = length(model.x);
    nu = length(model.u);

    if ~isempty(model.p)
        np = length(model.p);
    else
        np = 0;
    end

    if np > 0
        p0 = 1;
    end

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
    sim.code_gen_options.sens_forw_p = true;
    sim_solver = AcadosSimSolver(sim);

    sim_solver.set('x', x0);
    sim_solver.set('u', u);
    if np > 0
        sim_solver.set('p', p0);
    end

    sim_solver.solve();

    xn         = sim_solver.get('xn');
    S_forw_ind = sim_solver.get('S_forw');
    if np > 0
        S_p_ind = sim_solver.get('S_p');
    else
        S_p_ind = [];
    end

    if np > 0
        sim_solver.set('x', x0);
        sim_solver.set('u', u);
        sim_solver.set('p', p0);
        sim_solver.solve();
        xn_nom = sim_solver.get('xn');

        S_p_fd = zeros(nx, np);

        for jj = 1:np
            dp      = zeros(np, 1);
            dp(jj)  = 1.0;

            sim_solver.set('x', x0);
            sim_solver.set('u', u);
            sim_solver.set('p', p0 + FD_epsilon*dp);
            sim_solver.solve();

            xn_tmp = sim_solver.get('xn');
            S_p_fd(:, jj) = (xn_tmp - xn_nom) / FD_epsilon;
        end

        error_abs_Sp = max(max(abs(S_p_fd - S_p_ind)));
        disp(['error S_p (FD vs analytic):     ' num2str(error_abs_Sp)]);
        if error_abs_Sp > 1e-6
            error(['test_sens_forw_p FAIL: param sensitivities error too large for integrator ' method]);
        end
    else
        disp('Model has no parameters, skipping S_p test.');
    end
end

fprintf('\nTEST_PARAM_SENS (ERK+IRK): success!\n\n');