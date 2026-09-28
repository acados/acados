%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

%% test of native matlab interface (ERK/IRK + forward sens + param sens)

addpath('../pendulum_on_cart_model/');

for integrator = {'erk', 'irk'}

    %% integrator / method
    method = integrator{1};

    %% arguments
    compile_interface = 'auto';
    sens_forw   = 'true';
    sens_forw_p = 'true';   % param forward sensitivities
    jac_reuse   = 'true';
    num_stages  = 3;
    num_steps   = 4;
    newton_iter = 3;

    Ts         = 0.1;
    x0         = [1e-1; 1e0; 2e-1; 2e0];
    u          = 0;
    FD_epsilon = 1e-6;

    %% model
    model = pendulum_on_cart_model_with_param();
    model_name = ['pendulum_sens_p_' method];

    nx = model.nx;
    nu = model.nu;

    % detect parameters
    if isfield(model, 'sym_p')
        np = length(model.sym_p);
    else
        np = 0;
    end

    if np > 0
        p0 = 1;
    end

    %% acados sim model
    sim_model = acados_sim_model();
    sim_model.set('T', Ts);
    sim_model.set('name', model_name);

    sim_model.set('sym_x', model.sym_x);
    if isfield(model, 'sym_u')
        sim_model.set('sym_u', model.sym_u);
    end
    if isfield(model, 'sym_p')
        sim_model.set('sym_p', model.sym_p);
    end

    if (strcmp(method, 'erk'))
        sim_model.set('dyn_type', 'explicit');
        sim_model.set('dyn_expr_f', model.dyn_expr_f_expl);
    else
        sim_model.set('dyn_type', 'implicit');
        sim_model.set('dyn_expr_f', model.dyn_expr_f_impl);
        sim_model.set('sym_xdot', model.sym_xdot);
    end

    %% acados sim opts
    sim_opts = acados_sim_opts();
    sim_opts.set('compile_interface', compile_interface);
    sim_opts.set('num_stages', num_stages);
    sim_opts.set('num_steps', num_steps);
    sim_opts.set('newton_iter', newton_iter);
    sim_opts.set('method', method);
    sim_opts.set('sens_forw', sens_forw);
    sim_opts.set('sens_forw_p', sens_forw_p);
    sim_opts.set('jac_reuse', jac_reuse);

    %% acados sim
    sim_solver = acados_sim(sim_model, sim_opts);

    % set nominal state, input, parameter
    sim_solver.set('x', x0);
    sim_solver.set('u', u);
    if np > 0
        sim_solver.set('p', p0);
    end

    % solve once with analytic sensitivities on
    sim_solver.solve();

    xn         = sim_solver.get('xn');
    S_forw_ind = sim_solver.get('S_forw');
    if np > 0
        S_p_ind = sim_solver.get('S_p');
    else
        S_p_ind = [];
    end

    %% --- Param sensitivities S_p vs finite differences (p) ---

    if np > 0
        % Reset state, input, and parameter to nominal
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
        disp('Model has no parameters (np = 0), skipping S_p test.');
    end
end

fprintf('\nTEST_PARAM_SENS (ERK+IRK): success!\n\n');