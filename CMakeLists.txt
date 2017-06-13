#
#    This file is part of acados.
#
#    acados is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 3 of the License, or (at your option) any later version.
#
#    acados is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with acados; if not, write to the Free Software Foundation,
#    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#

cmake_minimum_required(VERSION 3.2)
project(acados LANGUAGES C CXX)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")
set(CMAKE_SKIP_RPATH TRUE) # Seems necessary for Mac OSX
set(CMAKE_VERBOSE_MAKEFILE ON)

message(STATUS "MAIN BEFORE EVERYTHING: ${CMAKE_CXX_FLAGS}")

# Set default compile flags for GCC
if(CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Werror -fPIC")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Werror -fPIC")
endif(CMAKE_COMPILER_IS_GNUCXX)

message(STATUS "MAIN AFTER CXX_FLAGS ASSIGNMENT: ${CMAKE_CXX_FLAGS}")

# Set up build type
if(NOT CMAKE_BUILD_TYPE)
    message("\n-- No build type chosen. Valid build types are:"
        "\n   Release (default), Debug, Test, NoExternalLibs")
    set(CMAKE_BUILD_TYPE Release)
endif()

message("\n\n --- Configuring acados with build type: ${CMAKE_BUILD_TYPE} ---\n\n")

if(CMAKE_BUILD_TYPE MATCHES Release OR CMAKE_BUILD_TYPE MATCHES NoExternalLibs)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -funroll-loops -fdiagnostics-show-option")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -funroll-loops -fdiagnostics-show-option")
endif()

message(STATUS "MAIN AFTER RELEASE: ${CMAKE_CXX_FLAGS}")

if(CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wno-unused-function -Wno-unused-parameter")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wno-unused-function -Wno-unused-parameter")
endif()

message(STATUS "MAIN AFTER DEBUG: ${CMAKE_CXX_FLAGS}")

if(CMAKE_BUILD_TYPE MATCHES Test)
    enable_testing(true)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

message(STATUS "MAIN AFTER TEST: ${CMAKE_CXX_FLAGS}")

# Set install and executable output directories
if(NOT ACADOS_INSTALL_DIR)
    set(ACADOS_INSTALL_DIR "$ENV{HOME}/local")
endif()
set(CMAKE_INSTALL_PREFIX ${ACADOS_INSTALL_DIR})
set(LIBRARY_OUTPUT_PATH "${CMAKE_INSTALL_PREFIX}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
message(STATUS "acados will be installed in ${CMAKE_INSTALL_PREFIX}")

message(STATUS "MAIN BEFORE SUBS: ${CMAKE_CXX_FLAGS}")

# Configure external libraries
set(EXTERNAL_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external)
add_subdirectory(${EXTERNAL_SRC_DIR})

message(STATUS "MAIN AFTER EXTERNAL: ${CMAKE_CXX_FLAGS}")

# Configure acados
set(ACADOS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/acados)
add_subdirectory(acados)

message(STATUS "MAIN AFTER ACADOS: ${CMAKE_CXX_FLAGS}")

# Add examples
set(EXAMPLES_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/examples)
add_subdirectory(${EXAMPLES_SRC_DIR})

message(STATUS "MAIN AFTER EXAMPLES: ${CMAKE_CXX_FLAGS}")

# Configure SWIG interface
if(NOT CMAKE_BUILD_TYPE MATCHES NoExternalLibs)
    add_subdirectory(swig)
endif()

# Generate test data and configure unit tests
if(CMAKE_BUILD_TYPE MATCHES Test)
    add_subdirectory(test)
endif()

# Cpplint
include(lint)
