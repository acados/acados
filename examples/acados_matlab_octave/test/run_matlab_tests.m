%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


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

%% run all tests
% IMPORTANT: the tests that are called here should NOT use clear all, but only call clear ocp_solver

test_names = [
    "test_code_reuse",
    "test_sim_code_reuse",
    "run_test_dim_check",
    "run_test_ocp_mass_spring",
    % "run_test_ocp_pendulum",
    "run_test_ocp_wtnx6",
    "run_test_sim_adj",
    "run_test_sim_dae",
    % "run_test_sim_forw",
    "run_test_sim_hess",
    "param_test",
    "test_conl_cost",
    "test_online_idxs_rev",
    "run_test_sim_sens_p",
    "test_slack_reformulation"
];

for k = 1:length(test_names)
    disp(strcat("running test ", test_names(k)));
    run(test_names(k))
end
