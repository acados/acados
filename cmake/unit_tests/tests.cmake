# Check if external libraries are present; these are needed for the unit tests
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/eigen/Eigen")
    message(FATAL_ERROR "The build type is ${CMAKE_BUILD_TYPE} (Test) but Eigen can not be found")
endif()
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/casadi-octave-v3.1.1")
    message(FATAL_ERROR "The build type is ${CMAKE_BUILD_TYPE} (Test) but casadi-octave-v3.1.1 can not be found")
endif()

# Automatic generation of unit testing data
include(unit_tests/test_sources)
add_custom_command(OUTPUT ${UNIT_TESTS_SRC_CASADI}
    COMMAND ${CMAKE_COMMAND} -P "${CMAKE_MODULE_PATH}/unit_tests/generate_test_data.cmake"
    COMMENT "Generating unit testing files")
add_custom_target(generate_test_data
    DEPENDS ${UNIT_TESTS_SRC_CASADI}
    COMMENT "Checking if regeneration of test data is needed")

# Unit test executable
add_executable(unit_tests ${UNIT_TESTS_SRC} ${UNIT_TESTS_SRC_SIM} ${UNIT_TESTS_SRC_OCP_QP} ${UNIT_TESTS_SRC_OCP_NLP} ${UNIT_TESTS_SRC_TEST_UTILS})
target_include_directories(unit_tests
    PRIVATE "${PROJECT_SOURCE_DIR}/external/eigen"
    )
target_link_libraries(unit_tests acados qpoases hpmpc blasfeo qpdunes m)
add_test(NAME unit_test_run COMMAND "${CMAKE_COMMAND}" -E chdir test ./unit_tests)
set_tests_properties(unit_test_run PROPERTIES DEPENDS generate_test_data)

if (EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
    find_package(FortranLibs REQUIRED)
    target_link_libraries(unit_tests ooqpgensparse ooqpsparse ooqpgondzio ooqpbase ma27 blas ${FORTRAN_LIBRARY} m)
endif ()

set_target_properties(unit_tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/test)

file(COPY "${PROJECT_SOURCE_DIR}/acados/sim/simplified/" DESTINATION "${PROJECT_BINARY_DIR}/test/simplified/")
