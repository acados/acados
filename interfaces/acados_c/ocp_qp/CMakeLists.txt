
set(INTERFACES_ACADOS_C_OCP_QP_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_full_condensing_solver.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_full_condensing.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_hpipm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_hpmpc.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_partial_condensing.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_qpdunes.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_sparse_solver.c)

target_sources(acados PRIVATE ${INTERFACES_ACADOS_C_OCP_QP_SRC})
