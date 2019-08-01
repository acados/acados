
set(INTERFACES_ACADOS_CPP_OCP_NLP_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/casadi_module.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/dynamic_loading.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_nlp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_nlp_solution.cpp)

target_sources(acados PRIVATE ${INTERFACES_ACADOS_CPP_OCP_NLP_SRC})
