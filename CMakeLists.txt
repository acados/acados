cmake_minimum_required(VERSION 2.6)
project(acados)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -v")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -Werror -fdiagnostics-show-option")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Werror")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -funroll-loops")
endif()

message("\n\n --- The CMake build type is ${CMAKE_BUILD_TYPE} ---\n\n")

if(CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG)
endif()

if(CMAKE_BUILD_TYPE MATCHES Profile)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -funroll-loops -fprofile-use")
endif()

if(CMAKE_BUILD_TYPE MATCHES Test)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
    enable_testing(true)
endif()

# Define include folders
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/external/")
include_directories("${PROJECT_SOURCE_DIR}/external/qpOASES/include/")
if(CMAKE_BUILD_TYPE MATCHES Test)
    if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/eigen/Eigen")
        message(FATAL_ERROR "The build type is ${CMAKE_BUILD_TYPE} (Test) but Eigen can not be found")
    endif()
    include_directories("${PROJECT_SOURCE_DIR}/external/eigen/")
endif()

# Define source files
include(acados_sources)
if(CMAKE_BUILD_TYPE MATCHES Test)
    include(generate_test_data)
    include(test_sources)
endif()

# Define external libraries
# BLASFEO
execute_process(COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/blasfeo
    make static_library RESULT_VARIABLE blasfeo_built)
if(NOT "${blasfeo_built}" STREQUAL "0")
    message(FATAL_ERROR "build of blasfeo failed with status ${blasfeo_built}!")
endif()
add_library(blasfeo STATIC IMPORTED)
set_property(TARGET blasfeo PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/external/blasfeo/libblasfeo.a")

# HPMPC
execute_process(COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/hpmpc
    make static_library RESULT_VARIABLE hpmpc_built)
if(NOT "${hpmpc_built}" STREQUAL "0")
    message(FATAL_ERROR "build of hpmpc failed with status ${hpmpc_built}!")
endif()
add_library(hpmpc STATIC IMPORTED)
set_property(TARGET hpmpc PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/external/hpmpc/libhpmpc.a")

# qpOASES
execute_process(COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/qpOASES
    make RESULT_VARIABLE qpoases_built)
if(NOT "${qpoases_built}" STREQUAL "0")
    message(FATAL_ERROR "build of qpoases failed with status ${qpoases_built}!")
endif()
add_library(qpoases STATIC IMPORTED)
set_property(TARGET qpoases PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/external/qpOASES/bin/libqpOASES_e.a")

# CasADi
if(BUILD_WITH_CASADI)
	execute_process(COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi mkdir build)
	execute_process(COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi cd build)
	execute_process(COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi/build cmake ..)
	execute_process(COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi/build
	    make RESULT_VARIABLE casadi_built)
	if(NOT "${casadi_built}" STREQUAL "0")
	    message(FATAL_ERROR "build of casadi failed with status ${casadi_built}!")
	endif()
endif()

# Define acados library
add_library(acados STATIC ${ACADOS_SRC} ${ACADOS_SRC_SIM} ${ACADOS_SRC_OCP_QP} ${ACADOS_SRC_OCP_NLP} ${ACADOS_SRC_UTILS})

# OpenBLAS installation directory
link_directories(/opt/openblas/lib)

# Define tests
add_executable(test_hpmpc ${TEST_HPMPC_SRC})
#target_link_libraries(test_hpmpc acados hpmpc blasfeo m)
target_link_libraries(test_hpmpc acados hpmpc blasfeo openblas m)

add_executable(test_condensing_qpoases ${TEST_CONDENSING_QPOASES_SRC})
#target_link_libraries(test_condensing_qpoases acados qpoases hpmpc blasfeo m)
target_link_libraries(test_condensing_qpoases acados qpoases hpmpc blasfeo openblas m)

add_executable(test_nmpc ${TEST_NMPC_SRC} ${TEST_CHEN_SRC})
#target_link_libraries(test_nmpc acados qpoases blasfeo m)
target_link_libraries(test_nmpc acados qpoases blasfeo openblas m)

add_executable(test_chain ${TEST_CHAIN_OCP_SRC} ${TEST_CHAIN_SRC} ${TEST_VDE_CHAIN_SRC} ${TEST_JAC_CHAIN_SRC})
#target_link_libraries(test_chain acados qpoases blasfeo m)
target_link_libraries(test_chain acados qpoases blasfeo openblas m)

if(CMAKE_BUILD_TYPE MATCHES Test)
    generate_test_data()
    add_executable(unit_tests ${UNIT_TESTS_SRC} ${UNIT_TESTS_SRC_OCP_QP} ${UNIT_TESTS_SRC_SIM} ${UNIT_TESTS_SRC_TEST_UTILS})
    target_link_libraries(unit_tests acados qpoases blasfeo hpmpc m)
	set_target_properties(unit_tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/test)
endif()

# Linter
find_package(PythonInterp)
if(PYTHONINTERP_FOUND)
    set(ACADOS_STYLE_FILTER
    	"-legal/copyright,"
    	"-readability/casting")
    set(FIND_COMMAND find acados test examples external/casadi/casadi/interfaces/acados -type f -name "*.c" ! \( -path "examples/casadi_chain/*" -o -path "test/sim/pendulum/casadi/*" \) -o -name "*.cpp" -o -name "*.h" -o -name "*.hpp")
    set(FIND_ACADOS_FILES ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${FIND_COMMAND})
    execute_process(COMMAND ${FIND_ACADOS_FILES} RESULT_VARIABLE RESULT OUTPUT_VARIABLE ACADOS_FILES_TO_LINT)
    string(REPLACE "\n" " " ACADOS_FILES_TO_LINT ${ACADOS_FILES_TO_LINT})
    separate_arguments(ACADOS_FILES_TO_LINT)
    set(LINT_COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${PYTHON_EXECUTABLE} ./cpplint.py --filter="${ACADOS_STYLE_FILTER}" --counting=detailed --extensions=c,cpp,h,hpp --linelength=100 ${ACADOS_FILES_TO_LINT})
    add_custom_target(lint ${LINT_COMMAND})
endif()
