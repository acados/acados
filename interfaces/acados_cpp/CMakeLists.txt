
set(INTERFACES_ACADOS_CPP_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_dimensions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp_solution.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp.cpp)

target_sources(acados PRIVATE ${INTERFACES_ACADOS_CPP_SRC})
