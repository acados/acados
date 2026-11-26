clear;
import casadi.*

angle = 1.0;
rotation = [cos(angle), -sin(angle); sin(angle), cos(angle)];
hessian = 1e12 * rotation * diag([1.0, 1e-8]) * rotation';
target = [2.0; 0.25];

model = AcadosModel();
model.name = 'hpipm_stats';
model.x = SX.sym('x', 2, 1);

ocp = AcadosOcp();
ocp.name = model.name;
ocp.model = model;
ocp.solver_options.N_horizon = 0;
ocp.solver_options.tf = 1.0;
ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.nlp_solver_max_iter = 10;
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.hpipm_mode = 'BALANCE';
ocp.solver_options.qp_solver_iter_max = 100;
ocp.solver_options.qp_solver_tol_stat = 1e-10;
ocp.solver_options.qp_solver_tol_eq = 1e-10;
ocp.solver_options.qp_solver_tol_ineq = 1e-10;
ocp.solver_options.qp_solver_tol_comp = 1e-10;
ocp.solver_options.nlp_solver_ext_qp_res = 1;
ocp.solver_options.hessian_approx = 'EXACT';
ocp.cost.cost_type_e = 'EXTERNAL';
ocp.model.cost_expr_ext_cost_e = 0.5 * (model.x - target)' * hessian * (model.x - target);
ocp.constraints.idxbx_e = [0; 1];
ocp.constraints.lbx_e = [-1.0; -1.0];
ocp.constraints.ubx_e = [1.0; 1.0];

ocp.simulink_opts = AcadosOcpSimulinkOptions();
input_names = properties(ocp.simulink_opts.inputs);
for i = 1:numel(input_names)
    ocp.simulink_opts.inputs.(input_names{i}) = 0;
end
output_names = properties(ocp.simulink_opts.outputs);
for i = 1:numel(output_names)
    ocp.simulink_opts.outputs.(output_names{i}) = 0;
end
ocp.simulink_opts.outputs.hpipm_lq_sum = 1;
ocp.simulink_opts.outputs.hpipm_itref_c_sum = 1;

ocp_solver = AcadosOcpSolver(ocp);
ocp_solver.solve();
if ocp_solver.get('status') ~= 0
    error('Reference SQP solve failed.');
end
expected_lq_sum = ocp_solver.get('hpipm_lq_sum');
expected_itref_c_sum = ocp_solver.get('hpipm_itref_c_sum');
if expected_lq_sum <= 0
    error('The BALANCE reference solve did not report an LQ fallback.');
end

original_directory = pwd;
cd(ocp.code_gen_options.code_export_directory);

simulink_model = [ocp.name '_ocp_solver_simulink_block'];
if bdIsLoaded(simulink_model)
    close_system(simulink_model, 0);
end
make_sfun;

simulink_model = [ocp.name '_ocp_solver_simulink_block'];
load_system(fullfile(ocp.code_gen_options.code_export_directory, [simulink_model '.slx']));
sfun_block = [simulink_model '/' ocp.name '_ocp_solver'];

add_block('simulink/Sinks/To Workspace', [simulink_model '/hpipm_lq_sum_sink'], ...
    'VariableName', 'hpipm_lq_sum_sim', 'SaveFormat', 'Timeseries');
add_block('simulink/Sinks/To Workspace', [simulink_model '/hpipm_itref_c_sum_sink'], ...
    'VariableName', 'hpipm_itref_c_sum_sim', 'SaveFormat', 'Timeseries');
add_line(simulink_model, [ocp.name '_ocp_solver/1'], 'hpipm_lq_sum_sink/1', 'autorouting', 'on');
add_line(simulink_model, [ocp.name '_ocp_solver/2'], 'hpipm_itref_c_sum_sink/1', 'autorouting', 'on');

sim_output = sim(simulink_model, 'StopTime', '1', 'ReturnWorkspaceOutputs', 'on');
lq_output = sim_output.get('hpipm_lq_sum_sim');
itref_c_output = sim_output.get('hpipm_itref_c_sum_sim');
fprintf('hpipm_lq_sum expected: %.16g, Simulink samples: %s\n', ...
    expected_lq_sum, mat2str(lq_output.Data(:)', 16));
fprintf('hpipm_itref_c_sum expected: %.16g, Simulink samples: %s\n', ...
    expected_itref_c_sum, mat2str(itref_c_output.Data(:)', 16));
if isempty(lq_output.Data) || abs(lq_output.Data(1) - expected_lq_sum) > 1e-10
    error('First Simulink hpipm_lq_sum sample does not match the reference solve.');
end
if isempty(itref_c_output.Data) || abs(itref_c_output.Data(1) - expected_itref_c_sum) > 1e-10
    error('First Simulink hpipm_itref_c_sum sample does not match the reference solve.');
end

close_system(simulink_model, 0);
cd(original_directory);
rmdir(ocp.code_gen_options.code_export_directory, 's');
disp('HPIPM Simulink output ports match the reference SQP solve.');
