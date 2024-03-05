%% Simulink example
%
clear all; clc;

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


%% Run the models
try
    sim('simulink_model_integrator.slx');
    cd ..
catch
    cd ..
    error('Simulink integrator example failed')
end

try
    cd c_generated_code
    sim('simulink_model_closed_loop.slx');
    cd ..
catch
    cd ..
    error('Simulink closed loop example failed')
end
