%% Simulink example
%

%% Run matlab example
%
example_closed_loop;


%% Render templated Code for the model contained in ocp object
%
ocp.generate_c_code;

%% Compile Sfunctions
cd c_generated_code

make_sfun; % ocp solver
make_sfun_sim; % integrator


%% Copy Simulink example blocks into c_generated_code
source_folder = fullfile(pwd, '..');
target_folder = pwd;
copyfile( fullfile(source_folder, 'simulink_model_closed_loop.slx'), target_folder );
copyfile( fullfile(source_folder, 'wind0_ref.mat'), target_folder );
copyfile( fullfile(source_folder, 'windN_ref.mat'), target_folder );
copyfile( fullfile(source_folder, 'y_ref.mat'), target_folder );
copyfile( fullfile(source_folder, 'y_e_ref.mat'), target_folder );


%% Open Simulink example blocks
open_system(fullfile(target_folder, 'simulink_model_closed_loop'))


%%
disp('Press play in Simulink!');
