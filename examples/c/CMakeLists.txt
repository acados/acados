# Define sources for models used in the examples
file(GLOB CHAIN_MODEL_SRC ${EXAMPLES_SRC_DIR}/c/chain_model/chain_model.c
    ${EXAMPLES_SRC_DIR}/c/chain_model/vde*.c
    ${EXAMPLES_SRC_DIR}/c/chain_model/jac*.c)

set(CHEN_MODEL_SRC ${EXAMPLES_SRC_DIR}/c/Chen_model/chen_model.c)
file(GLOB PENDULUM_MODEL_SRC ${EXAMPLES_SRC_DIR}/c/pendulum_model/pendulum_model.c
    ${EXAMPLES_SRC_DIR}/c/pendulum_model/vde_forw_pendulum.c
    ${EXAMPLES_SRC_DIR}/c/pendulum_model/jac_pendulum.c)

# Define examples
add_executable(mass_spring_hpmpc_example ${EXAMPLES_SRC_DIR}/c/mass_spring_hpmpc.c)
target_link_libraries(mass_spring_hpmpc_example acados)
add_test(mass_spring_hpmpc_example mass_spring_hpmpc_example)

# add_executable(mass_spring_partial_tightening_example
#    ${EXAMPLES_SRC_DIR}/c/mass_spring_partial_tightening.c
# )
# target_link_libraries(mass_spring_partial_tightening_example acados)
# add_test(
#     NAME mass_spring_partial_tightening_example
#     COMMAND "${CMAKE_COMMAND}" -E chdir ${EXECUTABLE_OUTPUT_PATH} ./mass_spring_partial_tightening_example
# )

#add_executable(nonlinear_chain_ocp_nlp # for sumulink interface
#  ${EXAMPLES_SRC_DIR}/c/nonlinear_chain_ocp_nlp.cpp
#  ${EXAMPLES_SRC_DIR}/c/chain_model/chain_model.c
#  ${EXAMPLES_SRC_DIR}/c/chain_model/jac_chain_nm2.c
#  ${EXAMPLES_SRC_DIR}/c/chain_model/jac_chain_nm3.c
#  ${EXAMPLES_SRC_DIR}/c/chain_model/jac_chain_nm4.c
#  ${EXAMPLES_SRC_DIR}/c/chain_model/jac_chain_nm5.c
#  ${EXAMPLES_SRC_DIR}/c/chain_model/jac_chain_nm6.c
#  ${EXAMPLES_SRC_DIR}/c/chain_model/jac_chain_nm7.c
#  ${EXAMPLES_SRC_DIR}/c/chain_model/jac_chain_nm8.c
#  ${EXAMPLES_SRC_DIR}/c/chain_model/jac_chain_nm9.c

#  ${EXAMPLES_SRC_DIR}/c/chain_model/vde_chain_nm2.c
#  ${EXAMPLES_SRC_DIR}/c/chain_model/jac_chain_nm2.c
#)
#target_link_libraries(nonlinear_chain_ocp_nlp acados)

add_executable(mass_spring_qpoases_example ${EXAMPLES_SRC_DIR}/c/mass_spring_qpoases.c)
target_link_libraries(mass_spring_qpoases_example acados)
add_test(mass_spring_qpoases_example mass_spring_qpoases_example)

add_executable(chen_nmpc_qpoases_example
    ${EXAMPLES_SRC_DIR}/c/chen_nmpc_qpoases.c ${CHEN_MODEL_SRC})
target_link_libraries(chen_nmpc_qpoases_example acados)
add_test(chen_nmpc_qpoases_example chen_nmpc_qpoases_example)

add_executable(nonlinear_chain_sqp_example
    ${EXAMPLES_SRC_DIR}/c/nonlinear_chain_sqp.c ${CHAIN_MODEL_SRC})
target_link_libraries(nonlinear_chain_sqp_example acados)
# Following test takes too long!
# add_test(
#    NAME nonlinear_chain_sqp_example
#    COMMAND "${CMAKE_COMMAND}" -E chdir ${EXECUTABLE_OUTPUT_PATH} ./nonlinear_chain_sqp_example
# )

# add_executable(pendulum_hpmpc_libstr_example ${EXAMPLES_SRC_DIR}/c/pendulum_hpmpc_libstr.c ${PENDULUM_MODEL_SRC})
# target_link_libraries(pendulum_hpmpc_libstr_example acados)
# add_test(NAME pendulum_hpmpc_libstr_example COMMAND "${CMAKE_COMMAND}" -E chdir ${EXECUTABLE_OUTPUT_PATH} ./pendulum_hpmpc_libstr_example)

add_executable(pendulum_partial_tightening_example
    ${EXAMPLES_SRC_DIR}/c/pendulum_partial_tightening.c ${PENDULUM_MODEL_SRC})
target_link_libraries(pendulum_partial_tightening_example acados)
add_test(pendulum_partial_tightening_example pendulum_partial_tightening_example)