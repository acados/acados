%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

clear VARIABLES

addpath('../pendulum_dae/');

i_method = 0;
for integrator = {'GNSF', 'IRK'}
    i_method = i_method + 1;
    method = integrator{1};

    sens_forw = true;
    jac_reuse = false;
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

    model = get_pendulum_dae_model();
    nx = length(model.x);
    nu = length(model.u);
    nz = length(model.z);

    sim = AcadosSim();
    sim.model = model;
    sim.model.name = model_name;
    sim.solver_options.Tsim = 0.1;
    sim.solver_options.integrator_type = method;
    sim.solver_options.num_stages = num_stages;
    sim.solver_options.num_steps = num_steps;
    sim.solver_options.newton_iter = newton_iter;
    sim.solver_options.sens_forw = sens_forw;
    sim.solver_options.sens_adj = true;
    sim.solver_options.sens_algebraic = true;
    sim.solver_options.output_z = true;
    sim.solver_options.jac_reuse = jac_reuse;
    sim_solver = AcadosSimSolver(sim);
    N_sim = 100;

    sim_solver.set('x', x0);
    sim_solver.set('u', u);

    x_sim = zeros(nx, N_sim+1);
    x_sim(:,1) = x0;

    tic
    for ii=1:N_sim

        sim_solver.set('x', x_sim(:,ii));
        sim_solver.set('u', u);
    
        sim_solver.set('seed_adj', ones(nx,1));
    
        if strcmp(method, 'IRK')
            sim_solver.set('xdot', zeros(nx,1));
            sim_solver.set('z', zeros(nz,1));
        elseif strcmp(method, 'GNSF')
            n_out = sim_solver.sim.model.gnsf_model.dims.nout;
            sim_solver.set('phi_guess', zeros(n_out,1));
        end

        sim_solver.solve();
    
        x_sim(:,ii+1) = sim_solver.get('xn');

    end
	S_forw = sim_solver.get('S_forw');
    S_adj = sim_solver.get('S_adj')';
    z = sim_solver.get('zn')';
    S_alg = sim_solver.get('S_algebraic');

    required_accuracy = 1e-13;
    if i_method == 1
        x_sim_ref = x_sim;
        S_forw_ref = S_forw;
        S_adj_ref = S_adj;
        z_ref = z;
        S_alg_ref = S_alg;
    else
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
