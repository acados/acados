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


%% create ocp solver
ocp = acados_ocp(ocp_model, ocp_opts, simulink_opts);

% solver initial guess
x_traj_init = rand(nx, N+1);
u_traj_init = zeros(nu, N);
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

status = ocp.get('status'); % 0 - success
ocp.print('stat')

%% simulink test
cd c_generated_code
make_sfun; % ocp solver
cd ..;

n_sim = 3;


%% Test Simulink example block
for itest = [1, 2, 3]
    if itest == 1
        % always reinitialize
        reset_value = 0;
        ignore_inits_value = 0;
    elseif itest == 2
        % always reset and initialize
        reset_value = 1;
        ignore_inits_value = 0;
    elseif itest == 3
        % always reset
        reset_value = 1;
        ignore_inits_value = 1;
    end
    out_sim = sim('initialization_test_simulink', 'SaveOutput', 'on');
    fprintf('\nSuccessfully ran simulink block with reset_value %d ignore_inits_value %d.\n\n', reset_value, ignore_inits_value);

    % Evaluation
    fprintf('\nTest results on SIMULINK simulation.\n')

    disp('checking KKT residual')
    kkt_signal = out_sim.logsout.getElement('KKT_residual');
    if any(kkt_signal.Values.data > 1e-6)
        disp('failed');
        quit(1);
    end

    sqp_iter_signal = out_sim.logsout.getElement('sqp_iter');
    sqp_iter_simulink = sqp_iter_signal.Values.Data;
    disp('checking SQP iter, QP should take 1 SQP iter.')
    if any(sqp_iter_simulink ~= 1)
        disp('failed');
        quit(1);
    end

    status_signal = out_sim.logsout.getElement('status');
    disp('checking status.')
    if any(status_signal.Values.Data)
        disp('failed. got status values:');
        disp(status_signal.Values.Data);
        % quit(1);
    end

    utraj_signal = out_sim.logsout.getElement('utraj');
    u_simulink = utraj_signal.Values.Data(1, :);
    disp('checking u values.')
    if any(abs(u_simulink(:) - utraj(:)) > 1e-8)
        disp('failed');
        quit(1);
    end

    xtraj_signal = out_sim.logsout.getElement('xtraj');
    x_simulink = xtraj_signal.Values.Data(1, :);
    disp('checking x values.')
    if any(abs(x_simulink(:) - xtraj(:)) > 1e-8)
        disp('failed');
        quit(1);
    end
end
%% Run with different initialization
% always reset and ignore initializations
reset_value = 1;
ignore_inits_value = 1;
out_sim = sim('initialization_test_simulink', 'SaveOutput', 'on');
disp('successfully ran simulink_model_advanced_closed_loop');

sqp_iter_signal = out_sim.logsout.getElement('sqp_iter');
sqp_iter_simulink = sqp_iter_signal.Values.Data;
disp('checking SQP iter, QP should take 1 SQP iter.')
if any(sqp_iter_simulink ~= 1)
    disp('failed');
    quit(1);
end

status_signal = out_sim.logsout.getElement('status');
disp('checking status.')
if any(status_signal.Values.Data)
    disp('failed. got status values:');
    disp(status_signal.Values.Data);
    quit(1);
end

%% Run with different initialization
% dont reset and ignore initializations -> only in first instance an SQP
% iteration is needed
reset_value = 0;
ignore_inits_value = 1;
out_sim = sim('initialization_test_simulink', 'SaveOutput', 'on');
disp('successfully ran simulink_model_advanced_closed_loop');

sqp_iter_signal = out_sim.logsout.getElement('sqp_iter');
sqp_iter_simulink = sqp_iter_signal.Values.Data;
disp('checking SQP iter, got')
disp(sqp_iter_simulink)

fprintf('should require one SQP iteration in first instance.\n')
if sqp_iter_simulink(1) ~= 1
    disp('failed');
    quit(1);
end
fprintf('should require 0 SQP iterations if initialized at solution.\n')
if any(sqp_iter_simulink(2:n_sim))
    disp('failed');
    quit(1);
end

status_signal = out_sim.logsout.getElement('status');
disp('checking status.')
if any(status_signal.Values.Data)
    disp('failed');
    quit(1);
end
