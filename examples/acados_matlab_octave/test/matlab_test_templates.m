%% run all Octave tests
test_names = ["run_test_template_ocp_disc_dyn",
"run_test_template_ocp_linear_dae",
"run_test_template_ocp_pendulum_auto",
"run_test_template_ocp_pendulum_exact_hess",
"run_test_template_ocp_pendulum_ext_cost",
"run_test_template_ocp_pendulum_gnsf",
% "run_test_template_ocp_pendulum_nls"
];

for k = 1:length(test_names)
    disp(strcat("running test ", test_names(k)));
    run(test_names(k))
end
