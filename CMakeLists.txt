cmake_minimum_required(VERSION 3.2)
project(acados LANGUAGES C CXX)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -Werror -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Werror -fPIC")

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type chosen. Valid build types are Release (default), Debug, Profile, Test")
    set(CMAKE_BUILD_TYPE Release)
endif()

message("\n\n --- Configuring acados with build type: ${CMAKE_BUILD_TYPE} ---\n\n")

if(CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -funroll-loops -fdiagnostics-show-option")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -funroll-loops -fdiagnostics-show-option")
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wno-unused-function -Wno-unused-parameter")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wno-unused-function -Wno-unused-parameter")
endif()

if(CMAKE_BUILD_TYPE MATCHES Profile)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -funroll-loops -fprofile-use")
	set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -O3 -funroll-loops -fprofile-use")
endif()

if(CMAKE_BUILD_TYPE MATCHES Test)
    enable_testing(true)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

# Set install and executable output directories
if (NOT ACADOS_INSTALL_DIR)
    set(ACADOS_INSTALL_DIR "$ENV{HOME}/local")
endif ()
set(CMAKE_INSTALL_PREFIX ${ACADOS_INSTALL_DIR})
set(LIBRARY_OUTPUT_PATH "${CMAKE_INSTALL_PREFIX}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
message(STATUS "acados will be installed in ${CMAKE_INSTALL_PREFIX}")

# Configure external libraries
set(EXTERNAL_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external)
add_subdirectory(${EXTERNAL_SRC_DIR})

# Configure acados
set(ACADOS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/acados)
add_subdirectory(acados)

# Add examples
set(EXAMPLES_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/examples)
add_subdirectory(${EXAMPLES_SRC_DIR})

# Configure SWIG interface
add_subdirectory(swig)

# Generate test data and configure unit tests
if(CMAKE_BUILD_TYPE MATCHES Test)
    include(unit_tests/tests)
endif()

# Cpplint
include(lint)
