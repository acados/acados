# Check if external libraries are present; these are needed for the unit tests
find_package(Eigen3)
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/casadi-octave-v3.1.1")
    message(FATAL_ERROR "The build type is ${CMAKE_BUILD_TYPE} (Test) but casadi-octave-v3.1.1 can not be found")
endif()

add_subdirectory(ocp_nlp)
add_subdirectory(ocp_qp)
add_subdirectory(sim)

set(TEST_OCP_NLP_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_nlp/chain/test_chain.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_nlp/chain/chain_model.c
)

set(TEST_OCP_QP_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp/test_condensing.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp/condensing_test_helper.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp/test_qpsolvers.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp/../test_utils/read_ocp_qp_in.c
)

set(TEST_SIM_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/sim/pendulum/test_pendulum.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sim/pendulum/casadi/casadi_pendulum.c
)

set(TEST_UTILS_SRC
    ${PROJECT_SOURCE_DIR}/test/test_utils/read_matrix.cpp
    ${PROJECT_SOURCE_DIR}/test/test_utils/zeros.cpp
) # TODO (dimitris): eliminate this at some point

# Unit test executable
add_executable(unit_tests
    ${PROJECT_SOURCE_DIR}/test/all_tests.cpp
    ${TEST_OCP_NLP_SRC}
    $<TARGET_OBJECTS:ocp_nlp_gen>
    ${TEST_OCP_QP_SRC}
    $<TARGET_OBJECTS:ocp_qp_gen>
    ${TEST_SIM_SRC}
    $<TARGET_OBJECTS:sim_gen>
    ${TEST_UTILS_SRC}
)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()
target_include_directories(unit_tests PRIVATE "${EXTERNAL_SRC_DIR}/eigen")
target_link_libraries(unit_tests acados)
set_target_properties(unit_tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/test)
if(EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
    target_compile_definitions(unit_tests PRIVATE OOQP)
endif()

# Add as test in ctest
add_test(NAME unit_tests COMMAND "${CMAKE_COMMAND}" -E chdir ${CMAKE_BINARY_DIR}/test ./unit_tests -a)

file(COPY "${PROJECT_SOURCE_DIR}/acados/sim/simplified/" DESTINATION "${PROJECT_BINARY_DIR}/test/simplified/")
