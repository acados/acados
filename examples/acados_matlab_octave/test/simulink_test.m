clear all
addpath(fullfile('..','getting_started'));

% extensive_example_ocp;


%% get available simulink_opts with default options
simulink_opts = get_acados_simulink_opts;

% manipulate simulink_opts

% inputs
simulink_opts.inputs.cost_W_0 = 1;
simulink_opts.inputs.cost_W = 1;
simulink_opts.inputs.cost_W_e = 1;
simulink_opts.inputs.x_init = 1;
simulink_opts.inputs.reset_solver = 1;

% outputs
simulink_opts.outputs.utraj = 1;
simulink_opts.outputs.xtraj = 1;

simulink_opts.samplingtime = '-1';
    % 't0' (default) - use time step between shooting node 0 and 1
    % '-1' - inherit sampling time from other parts of simulink model

% set time step for code generated integrator - default is length of first
% time step of ocp object

%% run minimal example
minimal_example_ocp;


%% Compile Sfunctions
cd c_generated_code

make_sfun_sim; % integrator
make_sfun; % ocp solver


%% Copy Simulink example blocks into c_generated_code
source_folder = fullfile(pwd, '..', '..', 'getting_started');
target_folder = pwd;
copyfile( fullfile(source_folder, 'simulink_model_advanced_closed_loop.slx'), target_folder );

%% Run Simulink example block
out_sim = sim('simulink_model_advanced_closed_loop', 'SaveOutput', 'on');

disp('successfully ran simulink_model_advanced_closed_loop');

cd ..

% simulink_u_traj_ref = out_sim.logsout{1}.Values.Data
% save('simulink_u_traj_ref.mat', 'simulink_u_traj_ref')
load('simulink_u_traj_ref.mat')

err_vs_ref = norm(simulink_u_traj_ref - out_sim.logsout{1}.Values.Data);
TOL = 1e-8;
disp(['Simulink: Norm of control traj. wrt. reference solution is: ',...
    num2str(err_vs_ref, '%e'), ' test TOL = ', num2str(TOL)]);

if err_vs_ref > TOL
    quit(1)
end



