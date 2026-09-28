%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


% NOTE: This example requires CasADi version 3.7 or later.
% Furthermore, this example requires additional flags for the CasADi code generation,
% cf. the solver option `ext_fun_compile_flags`

import casadi.*
lut = true;
use_p_global = true;
blazing = true;
fprintf('\n\nRunning Simulink test with lut=%d, use_p_global=%d, blazing=%d\n', lut, use_p_global, blazing);

% Create p_global parameters
[p_global, m, l, coefficients, ~, knots, p_global_values] = create_p_global(lut);

% OCP formulation
ocp = create_ocp_formulation_without_opts(p_global, m, l, coefficients, knots, lut, use_p_global, p_global_values, blazing);
ocp = set_solver_options(ocp);
ocp.model.name = ['sl_blz_' mat2str(blazing) '_pglbl_' mat2str(use_p_global) '_lut_' mat2str(lut)];
ocp.json_file = [ ocp.model.name '.json'];
% Simulink options
simulink_opts = get_acados_simulink_opts();  % should be AcadosOcpSimulinkOptions(), old function used to test backwards compatibility.
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
ocp.name = ocp.model.name;

% OCP solver
ocp_solver = AcadosOcpSolver(ocp);

%% MATLAB test solve
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
cd(ocp.code_gen_options.code_export_directory);
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
    error('xtraj values in SIMULINK and MATLAB should match.')
end
if norm(utraj_val - utraj) > 1e-8
    error('utraj values in SIMULINK and MATLAB should match.')
end

disp('Simulink p_global test: got matching trajectories in MATLAB and Simulink!')