%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

clear all;
check_acados_requirements()

import casadi.*

%%
N = 20; % number of discretization steps
nx = 3;
nu = 3;
np = 100;
[ocp, x0] = create_parametric_ocp_qp(N, np);

% NOTE: here we don't perform iterations and just test initialization
% functionality
ocp.solver_options.nlp_solver_max_iter = 0;
ocp.name = ocp.model.name;

% deactivate ports.
ocp.simulink_opts.inputs.lbx_0 = 0;
ocp.simulink_opts.inputs.ubx_0 = 0;
ocp.simulink_opts.inputs.lbx = 0;
ocp.simulink_opts.inputs.ubx = 0;
ocp.simulink_opts.inputs.reset_solver = 0;
ocp.simulink_opts.inputs.x_init = 0;
ocp.simulink_opts.inputs.u_init = 0;
ocp.simulink_opts.inputs.pi_init = 0;
ocp.simulink_opts.inputs.ignore_inits = 0;
ocp.simulink_opts.outputs.pi_all = 0;
ocp.simulink_opts.outputs.sqp_iter = 0;

% parameter ports
ocp.simulink_opts.inputs.parameter_traj = 1;
ocp.simulink_opts.outputs.parameter_traj = 1;

%% create ocp solver
ocp_solver = AcadosOcpSolver(ocp);

%% simulink test
cd(ocp.code_gen_options.code_export_directory);
make_sfun; % ocp solver
cd ..;

n_sim = 3;

p_values = 1:(np*(N+1));
out_sim = sim('parameter_test_simulink', 'SaveOutput', 'on');

status_signal = out_sim.logsout.getElement('status');
disp('checking status, should be 2 (max iter).')
if any(status_signal.Values.Data ~= 2)
    errror(['failed. got status values:' mat2str(status_signal.Values.Data)]);
end

element_names = out_sim.logsout.getElementNames();
parameter_traj_out_signal = out_sim.logsout.getElement('parameter_traj_out');
if any(parameter_traj_out_signal.Values.Data(1,:) ~= p_values)
    error('Setting parameters in Simulink does NOT work as expected.');
end
disp('Setting parameters in Simulink works as expected.');

