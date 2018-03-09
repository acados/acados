
set(INTERFACES_ACADOS_C_OCP_QP_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_full_condensing_solver_i.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_full_condensing_i.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_hpipm_i.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_partial_condensing_i.c
${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_sparse_solver_i.c)

if(ACADOS_WITH_QPDUNES)
    list(APPEND INTERFACES_ACADOS_C_OCP_QP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_qpdunes_i.c)
endif()

if(ACADOS_WITH_HPMPC)
    list(APPEND INTERFACES_ACADOS_C_OCP_QP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_hpmpc_i.c)
endif()

target_sources(acados PRIVATE ${INTERFACES_ACADOS_C_OCP_QP_SRC})
