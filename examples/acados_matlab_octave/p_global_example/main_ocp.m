%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


% NOTE: This example requires CasADi version 3.7 or later.
% Furthermore, this example requires additional flags for the CasADi code generation,
% cf. the solver option `ext_fun_compile_flags`


function main()

    import casadi.*

    %% Standard OCP compare blazing vs bspline, p global vs no p_global
    [state_trajectories_without_blazing_ref, t_tot_with_bspline_ref] = run_example_ocp(true, false, false);
    [state_trajectories_without_blazing, t_tot_with_bspline] = run_example_ocp(true, true, false);
    [state_trajectories_with_blazing_ref, t_tot_with_blazing_ref] = run_example_ocp(true, false, true);
    [state_trajectories_with_blazing, t_tot_with_blazing] = run_example_ocp(true, true, true);

    %% Timing comparison
    fprintf('\t\tbspline\t\tblazing\n');
    fprintf('ref\t\t%f \t%f\n', t_tot_with_bspline_ref, t_tot_with_blazing_ref);
    fprintf('p_global\t%f \t%f\n', t_tot_with_bspline, t_tot_with_blazing);

    %% Compare trajectories
    fprintf('max diff blazing with/without p_global %f\n', max(max(abs(state_trajectories_with_blazing_ref - state_trajectories_with_blazing))))
    fprintf('max diff blazing vs. bspline %f\n', max(max(abs(state_trajectories_with_blazing_ref - state_trajectories_without_blazing_ref))))

    %% Standard OCP without splines
    disp("Running OCP tests without splines.")
    [state_trajectories_no_lut_ref, ~] = run_example_ocp(false, false, true);
    [state_trajectories_no_lut, ~] = run_example_ocp(false, true, true);

    if ~(max(max(abs(state_trajectories_no_lut_ref - state_trajectories_no_lut))) < 1e-10)
        error("State trajectories with lut=false do not match.");
    end
end



function [state_trajectories, timing] = run_example_ocp(lut, use_p_global, blazing)

    import casadi.*

    fprintf('\n\nRunning example with lut=%d, use_p_global=%d, blazing=%d\n', lut, use_p_global, blazing);

    % Create p_global parameters
    [p_global, m, l, coefficients, ~, knots, p_global_values] = create_p_global(lut);

    % OCP formulation
    ocp = create_ocp_formulation_without_opts(p_global, m, l, coefficients, knots, lut, use_p_global, p_global_values, blazing);
    ocp = set_solver_options(ocp);
    ocp.model.name = ['ocp_blz_' mat2str(blazing) '_pglbl_' mat2str(use_p_global) '_lut_' mat2str(lut)];
    ocp.json_file = [ ocp.model.name '.json'];
    ocp.name = ocp.model.name;

    % OCP solver
    ocp_solver = AcadosOcpSolver(ocp);

    state_trajectories = [];  % only for testing purposes

    if use_p_global
        disp("Calling precompute.")
        tic
        ocp_solver.set_p_global_and_precompute_dependencies(p_global_values);
        toc
    end

    timing = 0;
    for i = 1:20
        ocp_solver.solve();
        state_trajectories = [state_trajectories; ocp_solver.get('x')];
        timing = timing + ocp_solver.get('time_lin');
    end

    % Plot results
    PLOT = false;

    if PLOT
        utraj = ocp_solver.get('u');
        xtraj = ocp_solver.get('x');
        plot_pendulum(ocp.solver_options.shooting_nodes, xtraj, utraj);
    end
end
