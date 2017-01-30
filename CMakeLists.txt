cmake_minimum_required(VERSION 2.8.11)
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
    set(CMAKE_C_FLAGS "-g -std=c99 -Wall -Wextra -Wno-unused-function")
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

# Process external libraries
include(external/openblas)
include(external/blasfeo)
include(external/hpmpc)
include(external/qpoases)
add_subdirectory(swig)

if (EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
    include(external/ma27)
    include(external/ooqp)
endif ()

# Define acados library
add_library(acados STATIC ${ACADOS_SRC})
add_dependencies(acados qpoases hpmpc blasfeo)

if (EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
    add_dependencies(ooqpbase ma27)
    add_dependencies(acados ooqpbase)
endif ()

# Define examples
add_executable(hpmpc_example ${HPMPC_EXAMPLE_SRC})
target_link_libraries(hpmpc_example acados hpmpc blasfeo openblas m)

add_executable(hpmpc_example_partial_tightening ${HPMPC_EXAMPLE_SRC})
target_link_libraries(hpmpc_example_partial_tightening acados hpmpc blasfeo openblas m)

add_executable(condensing_qpoases_example ${CONDENSING_QPOASES_EXAMPLE_SRC})
target_link_libraries(condensing_qpoases_example acados qpoases blasfeo hpmpc m)

add_executable(nmpc_example ${NMPC_EXAMPLE_SRC} ${CHEN_MODEL_SRC})
target_link_libraries(nmpc_example acados qpoases blasfeo m)

add_executable(chain_example ${CHAIN_EXAMPLE_SRC})
target_link_libraries(chain_example acados qpoases blasfeo openblas m)

add_executable(pendulum_example ${PENDULUM_EXAMPLE_SRC})
target_link_libraries(pendulum_example acados hpmpc blasfeo openblas m)

add_executable(pendulum_example_partial_tightening ${PENDULUM_EXAMPLE_PT_SRC})
target_link_libraries(pendulum_example_partial_tightening acados hpmpc blasfeo openblas m)

#add_executable(aircraft_example ${AIRCRAFT_EXAMPLE_SRC})
#target_link_libraries(aircraft_example acados hpmpc blasfeo openblas m)

# Automatic generation of unit testing data
if(CMAKE_BUILD_TYPE MATCHES Test)
    include(test_sources)
    add_custom_command(OUTPUT ${UNIT_TESTS_SRC_CASADI}
        COMMAND ${CMAKE_COMMAND} -P ../cmake/generate_test_data.cmake
        COMMENT "Generating unit testing files")
    add_custom_target(generate_test_data
        DEPENDS ${UNIT_TESTS_SRC_CASADI}
        COMMENT "Checking if regeneration of test data is needed")
    add_executable(unit_tests ${UNIT_TESTS_SRC} ${UNIT_TESTS_SRC_SIM} ${UNIT_TESTS_SRC_OCP_NLP} ${UNIT_TESTS_SRC_OCP_QP} ${UNIT_TESTS_SRC_TEST_UTILS})
    add_dependencies(unit_tests generate_test_data)

    if (EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
        find_library(FLIB gfortran HINTS /usr/lib/gcc/x86_64-linux-gnu/* /usr/local/lib/gcc/*)
        target_link_libraries(unit_tests acados qpoases blasfeo hpmpc ooqpgensparse ooqpsparse ooqpgondzio ooqpbase ma27 blas ${FLIB} m)
    else ()
        target_link_libraries(unit_tests acados qpoases blasfeo hpmpc m)
    endif ()
    set_target_properties(unit_tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/test)
endif()

# Linter
find_package(PythonInterp 3 REQUIRED)
set(FIND_FILES_TO_LINT find acados test examples swig -type f -name "*.c" -o -name "*.cpp" -o -name "*.h" -o -name "*.hpp" -o -name "*.i" )
set(FIND_FILES_TO_LINT ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${FIND_FILES_TO_LINT})
execute_process(COMMAND ${FIND_FILES_TO_LINT} OUTPUT_VARIABLE FILES_TO_LINT)
string(REPLACE "\n" " " FILES_TO_LINT ${FILES_TO_LINT})
separate_arguments(FILES_TO_LINT)
set(LINT_COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${PYTHON_EXECUTABLE} ./cpplint.py --counting=detailed --extensions=c,cpp,h,hpp,i ${FILES_TO_LINT})
add_custom_target(lint ${LINT_COMMAND})
