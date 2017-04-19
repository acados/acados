# Define sources for models used in the examples
file(GLOB CHAIN_MODEL_SRC ${ACADOS_EXAMPLES_DIR}/c/chain_model/chain_model.c
    ${ACADOS_EXAMPLES_DIR}/c/chain_model/vde*.c
    ${ACADOS_EXAMPLES_DIR}/c/chain_model/jac*.c
)
set(CHEN_MODEL_SRC ${ACADOS_EXAMPLES_DIR}/c/Chen_model/chen_model.c)
file(GLOB PENDULUM_MODEL_SRC ${ACADOS_EXAMPLES_DIR}/c/pendulum_model/pendulum_model.c
    ${ACADOS_EXAMPLES_DIR}/c/pendulum_model/vde_forw_pendulum.c
    ${ACADOS_EXAMPLES_DIR}/c/pendulum_model/jac_pendulum.c
)

# Define examples
add_executable(ocp_qp_hpmpc_example ${ACADOS_EXAMPLES_DIR}/c/ocp_qp_hpmpc.c)
target_link_libraries(ocp_qp_hpmpc_example acados hpmpc blasfeo openblas m)
add_test(NAME ocp_qp_hpmpc_example COMMAND "${CMAKE_COMMAND}" -E chdir ${EXECUTABLE_OUTPUT_PATH} ./ocp_qp_hpmpc_example)

# add_executable(ocp_qp_partial_tightening_example ${ACADOS_EXAMPLES_DIR}/c/ocp_qp_partial_tightening.c)
# target_link_libraries(ocp_qp_partial_tightening_example acados hpmpc blasfeo openblas m)
# add_test(NAME ocp_qp_partial_tightening_example COMMAND "${CMAKE_COMMAND}" -E chdir ${EXECUTABLE_OUTPUT_PATH} ./ocp_qp_partial_tightening_example)

add_executable(ocp_qp_qpoases_example ${ACADOS_EXAMPLES_DIR}/c/ocp_qp_qpoases.c)
target_link_libraries(ocp_qp_qpoases_example acados qpoases hpmpc blasfeo m)
add_test(NAME ocp_qp_qpoases_example COMMAND "${CMAKE_COMMAND}" -E chdir ${EXECUTABLE_OUTPUT_PATH} ./ocp_qp_qpoases_example)

add_executable(nmpc_qpoases_example ${ACADOS_EXAMPLES_DIR}/c/nmpc_qpoases.c ${CHEN_MODEL_SRC})
target_link_libraries(nmpc_qpoases_example acados qpoases blasfeo m)
add_test(NAME nmpc_qpoases_example COMMAND "${CMAKE_COMMAND}" -E chdir ${EXECUTABLE_OUTPUT_PATH} ./nmpc_qpoases_example)

add_executable(nonlinear_chain_sqp_example ${ACADOS_EXAMPLES_DIR}/c/nonlinear_chain_sqp.c ${CHAIN_MODEL_SRC})
target_link_libraries(nonlinear_chain_sqp_example acados qpoases blasfeo openblas m)
add_test(NAME nonlinear_chain_sqp_example COMMAND "${CMAKE_COMMAND}" -E chdir ${EXECUTABLE_OUTPUT_PATH} ./nonlinear_chain_sqp_example)

add_executable(pendulum_hpmpc_libstr_example ${ACADOS_EXAMPLES_DIR}/c/pendulum_hpmpc_libstr.c ${PENDULUM_MODEL_SRC})
target_link_libraries(pendulum_hpmpc_libstr_example acados hpmpc blasfeo openblas m)
add_test(NAME pendulum_hpmpc_libstr_example COMMAND "${CMAKE_COMMAND}" -E chdir ${EXECUTABLE_OUTPUT_PATH} ./pendulum_hpmpc_libstr_example)

add_executable(pendulum_partial_tightening_example ${ACADOS_EXAMPLES_DIR}/c/pendulum_partial_tightening.c ${PENDULUM_MODEL_SRC})
target_link_libraries(pendulum_partial_tightening_example acados hpmpc blasfeo openblas m)
add_test(NAME pendulum_partial_tightening_example COMMAND "${CMAKE_COMMAND}" -E chdir ${EXECUTABLE_OUTPUT_PATH} ./pendulum_partial_tightening_example)
