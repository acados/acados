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

    %% Multi-phase OCP without lut
    % these are commented for faster testing on CI.
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

    % [state_trajectories_with_lut_ref, ~, mocp_json] = run_example_mocp(true, false, true);
    % [state_trajectories_with_lut_ref, ~, mocp_json] = run_example_mocp(true, true, true);
    % [stat_traj_json_load_reuse, ~] = run_example_mocp_json_load(mocp_json, true);
    %
    [state_trajectories_with_lut_ref, ~, mocp_json] = run_example_mocp(true, false, true);
    [state_trajectories_with_lut, ~, mocp_json] = run_example_mocp(true, true, true);
    mocp_json = 'mocp_blz_true_pglbl_true_lut_true.json';
    [stat_traj_json_load, ~] = run_example_mocp_json_load(mocp_json, false);
    [stat_traj_json_load_reuse, ~] = run_example_mocp_json_load(mocp_json, true);

    if ~(max(max(abs(state_trajectories_with_lut - stat_traj_json_load))) < 1e-10)
        error("Results with loaded MOCP description does not match reference.");
    end
    if ~(max(max(abs(state_trajectories_with_lut_ref - state_trajectories_with_lut))) < 1e-10)
        error("State trajectories with lut=true do not match.");
    end
    if ~(max(max(abs(stat_traj_json_load_reuse - stat_traj_json_load))) < 1e-10)
        error("Results with loaded MOCP description does not match one with code reuse.");
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