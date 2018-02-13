
add_subdirectory(dense_qp)
add_subdirectory(ocp_qp)
add_subdirectory(sim)

set(INTERFACES_ACADOS_C_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/dense_qp_i.c
    ${CMAKE_CURRENT_SOURCE_DIR}/legacy_create_i.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_lin_i.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_i.c
    ${CMAKE_CURRENT_SOURCE_DIR}/options_i.c
${CMAKE_CURRENT_SOURCE_DIR}/sim_i.c)

target_sources(acados PRIVATE ${INTERFACES_ACADOS_C_SRC})
