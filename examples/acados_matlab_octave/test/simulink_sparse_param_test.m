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
np = 100;
[ocp_model, ocp_opts, simulink_opts, x0] = create_parametric_ocp_qp(N, np);
% NOTE: here we don't perform iterations and just test initialization
% functionality
ocp_opts.set('nlp_solver_max_iter', 0);

% deactivate ports.
simulink_opts.inputs.lbx_0 = 0;
simulink_opts.inputs.ubx_0 = 0;
simulink_opts.inputs.lbx = 0;
simulink_opts.inputs.ubx = 0;
simulink_opts.inputs.reset_solver = 0;
simulink_opts.inputs.x_init = 0;
simulink_opts.inputs.u_init = 0;
simulink_opts.inputs.pi_init = 0;
simulink_opts.inputs.ignore_inits = 0;
simulink_opts.outputs.pi_all = 0;
simulink_opts.outputs.sqp_iter = 0;

% parameter ports
simulink_opts.inputs.parameter_traj = 0;
simulink_opts.outputs.parameter_traj = 1;
% define multiple ports for sparse parameter update:
% Usage:
% add_sparse_param_port_simulink(simulink_opts, idx_p, port_name, stage_idx_0, stage_idx_e)
simulink_opts = add_sparse_param_port_simulink(simulink_opts, 0:7, 'first_8', 0, N);
simulink_opts = add_sparse_param_port_simulink(simulink_opts, [42, 43], 'p4243', 0, N);
simulink_opts = add_sparse_param_port_simulink(simulink_opts, 12, 'p12_stage3', 3, 3);
simulink_opts = add_sparse_param_port_simulink(simulink_opts, 12, 'p12_stage6', 6, 6);


%% create ocp solver
ocp_solver = acados_ocp(ocp_model, ocp_opts, simulink_opts);


%% test values
input_update_port_first8 = ones(8*(N+1) + 1, 1);
input_update_port_p4243 = [1; 8*ones(2*(N+1), 1)];
input_update_port_p12_stage3 = [1; 333];
input_update_port_p12_stage6 = [0; 5000];


ocp_solver.set_params_sparse(0:7, input_update_port_first8(2:9))
ocp_solver.set_params_sparse([42, 43], [8, 8]);
ocp_solver.set_params_sparse([12], input_update_port_p12_stage3(end), 3);



p_matlab = zeros(np, N+1);
for i=0:N
    p_matlab(:, i+1) = ocp_solver.get('p', i);
end

p_ref = zeros(np, N+1);
p_ref(1:8, :) = 1;
p_ref([43, 44], :) = 8;
p_ref([13], 3+1) = input_update_port_p12_stage3(end);


if any(any(p_ref ~= p_matlab))
    disp('Setting sparse parameters in Matlab does NOT work as expected.');
    quit(1);
else
    disp('Setting sparse parameters in Matlab works as expected.');
end

%% simulink test
cd c_generated_code
make_sfun; % ocp solver
cd ..;

n_sim = 3;


p_values = 1:(np*(N+1));
out_sim = sim('sparse_parameter_test_simulink', 'SaveOutput', 'on');

element_names = out_sim.logsout.getElementNames();
parameter_traj_out_signal = out_sim.logsout.getElement('parameter_traj_out');
p_simulink = reshape(parameter_traj_out_signal.Values.Data(1, :), np, N+1);

if any(any(p_simulink ~= p_matlab))
    disp('Setting sparse parameters in Simulink does NOT work as expected.');
    quit(1);
else
    disp('Setting sparse parameters in Simulink works as expected.');
end
