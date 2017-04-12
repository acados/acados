# Define examples
add_executable(hpmpc_example ${HPMPC_EXAMPLE_SRC})
target_link_libraries(hpmpc_example acados hpmpc blasfeo openblas m)

add_executable(hpmpc_example_partial_tightening ${HPMPC_EXAMPLE_SRC})
target_link_libraries(hpmpc_example_partial_tightening acados hpmpc blasfeo openblas m)

add_executable(condensing_qpoases_example ${CONDENSING_QPOASES_EXAMPLE_SRC})
target_link_libraries(condensing_qpoases_example acados qpoases hpmpc blasfeo m)

add_executable(nmpc_example ${NMPC_EXAMPLE_SRC} ${CHEN_MODEL_SRC})
target_link_libraries(nmpc_example acados qpoases blasfeo m)

add_executable(chain_example ${CHAIN_EXAMPLE_SRC})
target_link_libraries(chain_example acados qpoases blasfeo openblas m)

add_executable(pendulum_example ${PENDULUM_EXAMPLE_SRC})
target_link_libraries(pendulum_example acados hpmpc blasfeo openblas m)

add_executable(pendulum_example_partial_tightening ${PENDULUM_EXAMPLE_PT_SRC})
target_link_libraries(pendulum_example_partial_tightening acados hpmpc blasfeo openblas m)