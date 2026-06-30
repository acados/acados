
function [state_trajectories, timing, mocp_json] = run_example_mocp(lut, use_p_global, blazing)
    import casadi.*

    fprintf('\n\nRunning example with lut=%d, use_p_global=%d, blazing=%d\n', lut, use_p_global, blazing);

    % Create p_global parameters
    [p_global, m, l, coefficients, ~, knots, p_global_values] = create_p_global(lut);

    % MOCP formulation
    name = ['mocp_blz_' mat2str(blazing) '_pglbl_' mat2str(use_p_global) '_lut_' mat2str(lut)];
    mocp = create_mocp_formulation(p_global, m, l, coefficients, knots, lut, use_p_global, p_global_values, blazing, name);
    mocp.name = name;
    mocp.code_gen_options.json_file = [mocp.name '.json'];

    % MOCP solver
    mocp_solver = AcadosOcpSolver(mocp);

    state_trajectories = []; % only for testing purposes

    if use_p_global
        disp("Calling precompute.")
        tic
        mocp_solver.set_p_global_and_precompute_dependencies(p_global_values);
        toc
    end

    timing = 0;
    for i = 1:20
        mocp_solver.solve();
        state_trajectories = [state_trajectories; mocp_solver.get('x')];
        timing = timing + mocp_solver.get('time_lin');
    end

    % Plot results
    PLOT = false;

    if PLOT
        utraj = ocp_solver.get('u');
        xtraj = ocp_solver.get('x');
        plot_pendulum(ocp.solver_options.shooting_nodes, xtraj, utraj);
    end
    mocp_json = mocp.code_gen_options.json_file;
end


function mocp = create_mocp_formulation(p_global, m, l, coefficients, knots, lut, use_p_global, p_global_values, blazing, name)

    N_horizon_1 = 10;
    N_horizon_2 = 10;
    mocp = AcadosMultiphaseOcp([N_horizon_1, N_horizon_2]);
    ocp_phase_1 = create_ocp_formulation_without_opts(p_global, m, l, coefficients, knots, lut, use_p_global, p_global_values, blazing);
    ocp_phase_1.model.name = name;
    ocp_phase_2 = create_ocp_formulation_without_opts(p_global, m, l, coefficients, knots, lut, use_p_global, p_global_values, blazing);
    ocp_phase_2.model.name = name;

    mocp = set_solver_options(mocp);

    mocp.set_phase(ocp_phase_1, 1);
    mocp.set_phase(ocp_phase_2, 2);

    if use_p_global
        mocp.p_global_values = p_global_values;
    end
end
