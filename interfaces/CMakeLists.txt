
# add_test(NAME matlab_test_sim_pendulum 
#         COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/matlab_mex/pendulum_on_cart_model ${MATLAB_EXECUTABLE} -nodisplay -nodesktop -r  "run('./run_tests.m');")

add_test(NAME octave_test_sim_pendulum COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/matlab_mex/pendulum_on_cart_model octave-cli ./run_tests.m)