%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;


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


    %% Multi-phase OCP without lut
%     disp("Running MOCP tests without lut.")
%
%     [state_trajectories_no_lut_ref, ~] = run_example_mocp(false, false, true);
%     [state_trajectories_no_lut, ~] = run_example_mocp(false, true, true);
%
%     if ~(max(max(abs(state_trajectories_no_lut_ref - state_trajectories_no_lut))) < 1e-10)
%         error("State trajectories with lut=false do not match.");
%     end

    %% MOCP test with lut
    disp("Running MOCP tests with lut.")

    [state_trajectories_with_lut_ref, ~] = run_example_mocp(true, false, true);
    [state_trajectories_with_lut, ~] = run_example_mocp(true, true, true);

    if ~(max(max(abs(state_trajectories_with_lut_ref - state_trajectories_with_lut))) < 1e-10)
        error("State trajectories with lut=true do not match.");
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

function [state_trajectories, timing] = run_example_mocp(lut, use_p_global, blazing)
    import casadi.*

    fprintf('\n\nRunning example with lut=%d, use_p_global=%d, blazing=%d\n', lut, use_p_global, blazing);

    % Create p_global parameters
    [p_global, m, l, coefficients, ~, knots, p_global_values] = create_p_global(lut);

    % MOCP formulation
    name = ['mocp_blz_' mat2str(blazing) '_pglbl_' mat2str(use_p_global) '_lut_' mat2str(lut)];
    mocp = create_mocp_formulation(p_global, m, l, coefficients, knots, lut, use_p_global, p_global_values, blazing, name);
    mocp.name = name;
    mocp.json_file = [mocp.name '.json'];

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


function plot_pendulum(shooting_nodes, xtraj, utraj)
    figure; hold on;
    states = {'p', 'theta', 'v', 'dtheta'};
    for i=1:length(states)
        subplot(length(states), 1, i);
        plot(shooting_nodes, xtraj(i,:)); grid on;
        ylabel(states{i});
        xlabel('t [s]')
    end

    figure
    stairs(shooting_nodes, [utraj'; utraj(end)])
    ylabel('F [N]')
    xlabel('t [s]')
    grid on

    if is_octave()
        waitforbuttonpress;
    end
end