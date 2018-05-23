set(INTERFACES_ACADOS_C_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/external_function_interface.c
    ${CMAKE_CURRENT_SOURCE_DIR}/dense_qp_interface.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_nlp_interface.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_interface.c
    ${CMAKE_CURRENT_SOURCE_DIR}/condensing_interface.c
    ${CMAKE_CURRENT_SOURCE_DIR}/options_interface.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sim_interface.c)

target_sources(acados PRIVATE ${INTERFACES_ACADOS_C_SRC})