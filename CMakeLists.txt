cmake_minimum_required(VERSION 3.2)
project(acados)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -Werror -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Werror -fPIC")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message("\n\n --- The CMake build type is: ${CMAKE_BUILD_TYPE} ---\n\n")

if(CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -funroll-loops -fdiagnostics-show-option")
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wno-unused-function -Wno-unused-parameter")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wno-unused-function -Wno-unused-parameter")
endif()

if(CMAKE_BUILD_TYPE MATCHES Profile)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -funroll-loops -fprofile-use")
endif()

if(CMAKE_BUILD_TYPE MATCHES Test)
    enable_testing(true)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/local")
message(STATUS "Files will be installed in ${CMAKE_INSTALL_PREFIX}")

# Process external libraries
include(external/blasfeo)
include(external/hpmpc)
include(external/openblas)
include(external/qpdunes)
include(external/qpoases)
if(EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
    include(external/ma27)
    include(external/ooqp)
endif()

# Define acados library
include(acados_sources)
add_library(acados STATIC ${ACADOS_SRC})
target_include_directories(acados
    PUBLIC
        "${PROJECT_SOURCE_DIR}"
        "${PROJECT_SOURCE_DIR}/external/"
        "${PROJECT_SOURCE_DIR}/external/qpDUNES-dev/include/"
        "${PROJECT_SOURCE_DIR}/external/qpOASES/include/"
    )
add_dependencies(acados qpoases hpmpc blasfeo qpdunes)
if(EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
    target_compile_definitions(acados PRIVATE OOQP)
    add_dependencies(ooqpbase ma27)
    add_dependencies(acados ooqpbase)
endif()

include(examples)

add_subdirectory(swig)

# Generate test data and build unit tests
if(CMAKE_BUILD_TYPE MATCHES Test)
    include(unit_tests/tests)
endif()

include(lint)
