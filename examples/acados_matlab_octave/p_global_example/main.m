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


% NOTE: This example requires CasADi version nightly-se2 or later,
% as well as an installation of simde.
% Furthermore, this example requires additional flags for the CasADi code generation,
% cf. the solver option `ext_fun_compile_flags`

function main()

    import casadi.*
    % Standard OCP
    state_trajectories_no_lut_ref = run_example_ocp(false, false);
    state_trajectories_no_lut = run_example_ocp(false, true);

    if ~all(abs(state_trajectories_no_lut_ref - state_trajectories_no_lut) < 1e-10)
        error("State trajectories with lut=false do not match.");
    end

    state_trajectories_with_lut_ref = run_example_ocp(true, false);
    state_trajectories_with_lut = run_example_ocp(true, true);

    if ~all(abs(state_trajectories_with_lut_ref - state_trajectories_with_lut) < 1e-10)
        error("State trajectories with lut=true do not match.");
    end

    % Multi-phase OCP
    state_trajectories_no_lut_ref = run_example_mocp(false, false);
    state_trajectories_no_lut = run_example_mocp(false, true);

    if ~all(abs(state_trajectories_no_lut_ref - state_trajectories_no_lut) < 1e-10)
        error("State trajectories with lut=false do not match.");
    end

    state_trajectories_with_lut_ref = run_example_mocp(true, false);
    state_trajectories_with_lut = run_example_mocp(true, true);

    if ~all(abs(state_trajectories_with_lut_ref - state_trajectories_with_lut) < 1e-10)
        error("State trajectories with lut=true do not match.");
    end
    %% Simulink test
    if ~is_octave()
        run_example_ocp_simulink_p_global();
    end

end



function run_example_ocp_simulink_p_global()

    import casadi.*
    lut = true;
    use_p_global = true;
    fprintf('\n\nRunning example with lut=%d, use_p_global=%d\n', lut, use_p_global);

    % Create p_global parameters
    [p_global, m, l, C, p_global_values] = create_p_global(lut);

    % OCP formulation
    ocp = create_ocp_formulation_without_opts(p_global, m, l, C, lut, use_p_global, p_global_values);
    ocp = set_solver_options(ocp);

    % Simulink options
    simulink_opts = get_acados_simulink_opts();
    simulink_opts.inputs.p_global = 1;
    possible_inputs = fieldnames(simulink_opts.inputs);
    for i = 1:length(possible_inputs)
        simulink_opts.inputs.(possible_inputs{i}) = 0;
    end
    simulink_opts.inputs.lbx_0 = 1;
    simulink_opts.inputs.ubx_0 = 1;
    simulink_opts.inputs.p_global = 1;

    simulink_opts.outputs.xtraj = 1;
    simulink_opts.outputs.utraj = 1;
    simulink_opts.outputs.u0 = 0;
    simulink_opts.outputs.x1 = 0;

    ocp.simulink_opts = simulink_opts;

    % OCP solver
    ocp_solver = AcadosOcpSolver(ocp);

    %% Matlab test solve
    % test with ones such that update is necessary
    p_global_values_test = ones(size(p_global_values));
    if use_p_global
        ocp_solver.set_p_global_and_precompute_dependencies(p_global_values_test);
    end

    ocp_solver.solve();
    xtraj = ocp_solver.get('x');
    xtraj = xtraj(:)';
    utraj = ocp_solver.get('u');
    utraj = utraj(:)';


    %% build s funtion
    cd c_generated_code;
    make_sfun;
    cd ..;

    %% run simulink block
    out_sim = sim('p_global_simulink_test_block', 'SaveOutput', 'on');
    fprintf('\nSuccessfully ran simulink block');

    %% Evaluation
    fprintf('\nTest results on SIMULINK simulation.\n')

    disp('checking KKT residual')
    % kkt_signal = out_sim.logsout.getElement('KKT_residual');
    xtraj_signal = out_sim.logsout.getElement('xtraj');
    xtraj_val = xtraj_signal.Values.Data(1, :);
    utraj_signal = out_sim.logsout.getElement('utraj');
    utraj_val = utraj_signal.Values.Data(1, :);
    if norm(xtraj_val - xtraj) > 1e-8
        disp('error: xtraj values in SIMULINK and MATLAB should match.')
        quit(1);
    end
    if norm(utraj_val - utraj) > 1e-8
        disp('error: utraj values in SIMULINK and MATLAB should match.')
        quit(1);
    end
    disp('Simulink p_global test: got matching trajectories in Matlab and Simulink!')
end


function state_trajectories = run_example_ocp(lut, use_p_global)

    import casadi.*

    fprintf('\n\nRunning example with lut=%d, use_p_global=%d\n', lut, use_p_global);

    % Create p_global parameters
    [p_global, m, l, C, p_global_values] = create_p_global(lut);

    % OCP formulation
    ocp = create_ocp_formulation_without_opts(p_global, m, l, C, lut, use_p_global, p_global_values);
    ocp = set_solver_options(ocp);

    % OCP solver
    ocp_solver = AcadosOcpSolver(ocp);

    state_trajectories = [];  % only for testing purposes

    if use_p_global
        ocp_solver.set_p_global_and_precompute_dependencies(p_global_values);
    end

    for i = 1:20
        ocp_solver.solve();
        state_trajectories = [state_trajectories; ocp_solver.get('x')];
    end

    % Plot results
    PLOT = false;

    if PLOT
        utraj = ocp_solver.get('u');
        xtraj = ocp_solver.get('x');
        plot_pendulum(ocp.solver_options.shooting_nodes, xtraj, utraj);
    end
end

function state_trajectories = run_example_mocp(lut, use_p_global)
    import casadi.*

    fprintf('\n\nRunning example with lut=%d, use_p_global=%d\n', lut, use_p_global);

    % Create p_global parameters
    [p_global, m, l, C, p_global_values] = create_p_global(lut);

    % OCP formulation
    mocp = create_mocp_formulation(p_global, m, l, C, lut, use_p_global, p_global_values);

    % OCP solver
    mocp_solver = AcadosOcpSolver(mocp);

    state_trajectories = []; % only for testing purposes

    if use_p_global
        mocp_solver.set_p_global_and_precompute_dependencies(p_global_values);
    end

    for i = 1:20
        mocp_solver.solve();
        state_trajectories = [state_trajectories; mocp_solver.get('x')];
    end

    % Plot results
    PLOT = false;

    if PLOT
        utraj = ocp_solver.get('u');
        xtraj = ocp_solver.get('x');
        plot_pendulum(ocp.solver_options.shooting_nodes, xtraj, utraj);
    end
end


function ocp = set_solver_options(ocp)
    % set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
    ocp.solver_options.integrator_type = 'ERK';
    ocp.solver_options.print_level = 0;
    ocp.solver_options.nlp_solver_type = 'SQP_RTI';

    % set prediction horizon
    Tf = 1.0;
    N_horizon = 20;
    ocp.solver_options.tf = Tf;
    ocp.solver_options.N_horizon = N_horizon;

    % partial condensing
    ocp.solver_options.qp_solver_cond_N = 5;
    ocp.solver_options.qp_solver_cond_block_size = [3, 3, 3, 3, 7, 1];

    % NOTE: these additional flags are required for code generation of CasADi functions using casadi.blazing_spline
    ocp.solver_options.ext_fun_compile_flags = ['-I' casadi.GlobalOptions.getCasadiIncludePath ' -ffast-math -march=native '];
end

function mocp = create_mocp_formulation(p_global, m, l, C, lut, use_p_global, p_global_values)

    N_horizon_1 = 10;
    N_horizon_2 = 10;
    mocp = AcadosMultiphaseOcp([N_horizon_1, N_horizon_2]);
    ocp_phase_1 = create_ocp_formulation_without_opts(p_global, m, l, C, lut, use_p_global, p_global_values);
    ocp_phase_2 = create_ocp_formulation_without_opts(p_global, m, l, C, lut, use_p_global, p_global_values);

    mocp = set_solver_options(mocp);

    mocp.set_phase(ocp_phase_1, 1);
    mocp.set_phase(ocp_phase_2, 2);

    if use_p_global
        mocp.p_global_values = p_global_values;
    end
end


function [p_global, m, l, C, p_global_values] = create_p_global(lut)

    import casadi.*
    m = MX.sym('m');
    l = MX.sym('l');
    p_global = {m, l};
    p_global_values = [0.1; 0.8];

    if lut
        data = rand(7, 5); % Example data, replace with actual data
        C = MX.sym('C', numel(data), 1);
        p_global{end+1} = C;
        p_global_values = [p_global_values; data(:)];
    else
        C = [];
    end

    p_global = vertcat(p_global{:});
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