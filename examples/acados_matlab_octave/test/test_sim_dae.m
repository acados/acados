%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

%% test of native matlab interface
clear VARIABLES

addpath('../pendulum_dae/');

i_method = 0;
for integrator = {'irk_gnsf', 'irk'}
    i_method = i_method + 1;
    method = integrator{1};

    %% arguments
    compile_interface = 'auto'; % true, false
    gnsf_detect_struct = 'true'; % true, false
    % method = 'irk'; % irk, irk_gnsf, [erk]
    sens_forw = 'true'; % true, false
    jac_reuse = 'false'; % true, false
    num_stages = 3;
    num_steps = 3;
    newton_iter = 5;
    model_name = ['pend_dae_' method];

    length_pendulum = 5;
    alpha0 = .01;
    xp0 = length_pendulum * sin(alpha0);
    yp0 = - length_pendulum * cos(alpha0);
    x0 = [ xp0; yp0; alpha0; 0; 0; 0];

    u = 3.5;

    %% model
    model = pendulum_dae_model();
    % disp('state')
    % disp(model.sym_x)

    nx = length(model.sym_x);
    nu = length(model.sym_u);
    nz = length(model.sym_z);

    %% acados sim model
    sim_model = acados_sim_model();
    sim_model.set('name', model_name);
    sim_model.set('T', 0.1); % simulation time

    sim_model.set('sym_x', model.sym_x);
    if isfield(model, 'sym_u')
        sim_model.set('sym_u', model.sym_u);
    end
    if isfield(model, 'sym_p')
        sim_model.set('sym_p', model.sym_p);
    end

    % Note: DAEs can only be used with implicit integrator
    sim_model.set('dyn_type', 'implicit');
    sim_model.set('dyn_expr_f', model.expr_f_impl);
    sim_model.set('sym_xdot', model.sym_xdot);
    if isfield(model, 'sym_z')
        sim_model.set('sym_z', model.sym_z);
    end

    %% acados sim opts
    sim_opts = acados_sim_opts();
    sim_opts.set('compile_interface', compile_interface);
    sim_opts.set('num_stages', num_stages);
    sim_opts.set('num_steps', num_steps);
    sim_opts.set('newton_iter', newton_iter);
    sim_opts.set('method', method);
    sim_opts.set('sens_forw', sens_forw);
    sim_opts.set('sens_adj', 'true');
    sim_opts.set('sens_algebraic', 'true');
    sim_opts.set('output_z', 'true');
    sim_opts.set('sens_hess', 'false');
    sim_opts.set('jac_reuse', jac_reuse);
    if (strcmp(method, 'irk_gnsf'))
        sim_opts.set('gnsf_detect_struct', gnsf_detect_struct);
    end


    %% acados sim
    % create integrator
    sim_solver = acados_sim(sim_model, sim_opts);
    N_sim = 100;

    % set initial state
    sim_solver.set('x', x0);
    sim_solver.set('u', u);

    x_sim = zeros(nx, N_sim+1);
    x_sim(:,1) = x0;

    tic
    for ii=1:N_sim

        % set initial state
        sim_solver.set('x', x_sim(:,ii));
        sim_solver.set('u', u);
    
        % set adjoint seed
        sim_solver.set('seed_adj', ones(nx,1));
    
        % initialize implicit integrator
        if (strcmp(method, 'irk'))
            sim_solver.set('xdot', zeros(nx,1));
            sim_solver.set('z', zeros(nz,1));
        elseif (strcmp(method, 'irk_gnsf'))
            n_out = sim_solver.sim.model.gnsf_model.dims.nout;
            sim_solver.set('phi_guess', zeros(n_out,1));
        end

        % solve
        sim_solver.solve();
    
        % get simulated state
        x_sim(:,ii+1) = sim_solver.get('xn');

    end
	S_forw = sim_solver.get('S_forw');
    S_adj = sim_solver.get('S_adj')';
    z = sim_solver.get('zn')'; % approximate value of algebraic variables at start of simulation
    S_alg = sim_solver.get('S_algebraic'); % sensitivities of algebraic variables z

    required_accuracy = 1e-13;
    if i_method == 1
        % store solution as reference
        x_sim_ref = x_sim;
        S_forw_ref = S_forw;
        S_adj_ref = S_adj;
        z_ref = z;
        S_alg_ref = S_alg;
    else
%         S_alg
%         S_alg_ref

        err_x = norm(x_sim - x_sim_ref);
        err_S_forw = norm(S_forw - S_forw_ref);
        err_S_adj = norm(S_adj - S_adj_ref);
        err_z = norm(z - z_ref);
        err_S_alg = norm(S_alg - S_alg_ref);
        err = max([err_x, err_S_forw, err_S_adj, err_z, err_S_alg]);
        fprintf(['\nerr_x\t\t' num2str(err_x, '%e')]);
        fprintf(['\nerr_S_forw\t' num2str(err_S_forw, '%e')]);
        fprintf(['\nerr_S_adj\t' num2str(err_S_adj, '%e')]);
        fprintf(['\nerr_z\t\t' num2str(err_z, '%e')]);
        fprintf(['\nerr_S_alg\t' num2str(err_S_alg, '%e')]);

        if max(err > required_accuracy )
            error(strcat('test_sim_dae FAIL: error larger than required accuracy:',...
                num2str(required_accuracy), ' for integrator: ', method));
        end
    end


end

fprintf('\n\nTEST_SIM_DAE: success!\n\n');
