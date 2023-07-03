%% Simulink example
%
clear all;

%% Run minimal example
%
minimal_example_ocp;


%% Compile Sfunctions
cd c_generated_code

make_sfun; % ocp solver
make_sfun_sim; % integrator


%% Copy Simulink example blocks into c_generated_code
source_folder = fullfile(pwd, '..');
target_folder = pwd;
copyfile( fullfile(source_folder, 'simulink_model_integrator.slx'), target_folder );
copyfile( fullfile(source_folder, 'simulink_model_closed_loop.slx'), target_folder );


%% Open Simulink example blocks
open_system(fullfile(target_folder, 'simulink_model_integrator'))
open_system(fullfile(target_folder, 'simulink_model_closed_loop'))


%%
disp('Press play in Simulink!');
