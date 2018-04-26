
target_sources(acados PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/options.cpp)

add_subdirectory(ocp_nlp)
add_subdirectory(ocp_qp)
