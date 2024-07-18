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

%

clear all;
check_acados_requirements()

import casadi.*

%%
N = 20; % number of discretization steps
nx = 3;
nu = 3;
[ocp_model, ocp_opts, simulink_opts, x0] = create_ocp_qp_solver_formulation(N);
% NOTE: here we don't perform iterations and just test initialization
% functionality
ocp_opts.set('nlp_solver_max_iter', 0);


%% create ocp solver
ocp = acados_ocp(ocp_model, ocp_opts, simulink_opts);

% solver initial guess
x_traj_init = rand(nx, N+1);
u_traj_init = rand(nu, N);
pi_init = rand(nx, N);

%% call ocp solver
% update initial state
ocp.set('constr_x0', x0);

% set trajectory initialization
ocp.set('init_x', x_traj_init); % states
ocp.set('init_u', u_traj_init); % inputs
ocp.set('init_pi', pi_init); % multipliers for dynamics equality constraints

% solve
ocp.solve();
% get solution
utraj = ocp.get('u');
xtraj = ocp.get('x');
pi_all = ocp.get('pi');

if norm(pi_init - pi_all) > 1e-10
    disp('pi initialization in MEX failed')
end
if norm(utraj - u_traj_init) > 1e-10
    disp('u initialization in MEX failed')
end
if norm(xtraj - x_traj_init) > 1e-10
    disp('x initialization in MEX failed')
end

status = ocp.get('status'); % 0 - success
ocp.print('stat')

%% simulink test
cd c_generated_code
make_sfun; % ocp solver
cd ..;
n_sim = 3;

%% Test Simulink example block
for itest = [1, 2, 3, 4]
    if itest == 1
        % always reinitialize
        reset_value = 0;
        ignore_inits_value = 0;
    elseif itest == 2
        % always reset
        reset_value = 1;
        ignore_inits_value = 0;
    elseif itest == 3
        % dont reset and dont initialize
        reset_value = 0;
        ignore_inits_value = 1;
    elseif itest == 4
        % always reset and initialize
        reset_value = 1;
        ignore_inits_value = 1;
    end

    if (itest == 1 || itest == 2)
        u_expected = u_traj_init;
        x_expected = x_traj_init;
        pi_expected = pi_init;
    elseif (itest == 3)
        u_expected = 0 * u_traj_init;
        % Note: first solver call is initialized with x0 for whole horizon.
        x_expected = repmat(x0, 1, N+1);
        pi_expected = 0 * pi_init;
    elseif (itest == 4)
        u_expected = 0 * u_traj_init;
        pi_expected = 0 * pi_init;
        x_expected = 0 * x_traj_init;
    end

    out_sim = sim('initialization_test_simulink', 'SaveOutput', 'on');
    fprintf('\nSuccessfully ran simulink block with reset_value %d ignore_inits_value %d.\n\n', reset_value, ignore_inits_value);

    % Evaluation
    fprintf('\nTest results on SIMULINK simulation.\n')

    disp('checking KKT residual, should be zero just because not evaluated')
    kkt_signal = out_sim.logsout.getElement('KKT_residual');
    if any(kkt_signal.Values.data > 1e-6)
        disp('failed');
        quit(1);
    end

    sqp_iter_signal = out_sim.logsout.getElement('sqp_iter');
    disp('checking SQP iter, we set max iter to 0, so expect 0.')
    if any(sqp_iter_signal.Values.Data ~= 0)
        disp('failed');
        quit(1);
    end

    status_signal = out_sim.logsout.getElement('status');
    disp('checking status, should be 2 (max iter).')
    if any(status_signal.Values.Data ~= 2)
        disp('failed. got status values:');
        disp(status_signal.Values.Data);
        quit(1);
    end

    utraj_signal = out_sim.logsout.getElement('utraj');
    u_simulink = utraj_signal.Values.Data(1, :);
    disp('checking u values.')
    if any(abs(u_simulink(:) - u_expected(:)) > 1e-8)
        disp('failed');
        quit(1);
    end

    xtraj_signal = out_sim.logsout.getElement('xtraj');
    x_simulink = xtraj_signal.Values.Data(1, :);
    disp('checking x values, should match initialization.')
    if any(abs(x_simulink(:) - x_expected(:)) > 1e-8)
        disp('failed');
        quit(1);
    end

    pi_signal = out_sim.logsout.getElement('pi_all');
    pi_simulink = pi_signal.Values.Data(1, :);
    disp('checking pi values, should match initialization.')
    if any(abs(pi_simulink(:) - pi_expected(:)) > 1e-8)
        disp('failed');
        quit(1);
    end
end