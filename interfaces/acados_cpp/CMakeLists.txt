
target_sources(acados PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/code_generator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/function_generation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_bounds.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_dimensions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/integrator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/options.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/utils.cpp)

add_subdirectory(ocp_nlp)
add_subdirectory(ocp_qp)
# todo add sim

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp")
