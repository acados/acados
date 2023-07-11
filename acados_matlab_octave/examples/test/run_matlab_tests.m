assert(1+1==2)

disp('assertation works')

disp('checking environment variables')

disp('MATLABPATH')
disp(getenv('MATLABPATH'))

disp('MODEL_FOLDER')
disp(getenv('MODEL_FOLDER'))


disp('ENV_RUN')
disp(getenv('ENV_RUN'))

disp('LD_LIBRARY_PATH')
disp(getenv('LD_LIBRARY_PATH'))

disp('pwd')
disp(pwd)

disp('running tests')

%% run all Octave tests
test_names = ["run_test_dim_check",
"run_test_ocp_mass_spring",
% "run_test_ocp_pendulum",
"run_test_ocp_pendulum_dae",
"run_test_ocp_simple_dae",
"run_test_ocp_wtnx6",
% "run_test_sim_adj",
"run_test_sim_dae",
% "run_test_sim_forw",
"run_test_sim_hess",
"run_test_mhe_lorentz"
];

for k = 1:length(test_names)
    disp(strcat("running test ", test_names(k)));
    run(test_names(k))
end
