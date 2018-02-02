
add_subdirectory(dense_qp)
add_subdirectory(ocp_qp)
add_subdirectory(sim)

set(INTERFACES_ACADOS_C_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/dense_qp.c
    ${CMAKE_CURRENT_SOURCE_DIR}/legacy_create.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_lin.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp.c
    ${CMAKE_CURRENT_SOURCE_DIR}/options.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sim.c)

target_sources(acados PRIVATE ${INTERFACES_ACADOS_C_SRC})
