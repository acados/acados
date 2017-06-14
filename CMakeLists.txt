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
set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

# Set up build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release
        CACHE STRING "No build type chosen. Valid build types are: \
        Release (default), None, Debug, NoExternalLibs, RelWithDebInfo, MinSizeRel." FORCE)
endif()
message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")

enable_testing(true)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wno-unused-function -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wno-unused-function -Wno-unused-parameter")

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g")

# Set install directories
if(NOT ACADOS_INSTALL_DIR)
    set(ACADOS_INSTALL_DIR "$ENV{HOME}/local")
endif()
set(CMAKE_INSTALL_PREFIX ${ACADOS_INSTALL_DIR}
    CACHE STRING "Install path prefix, prepended onto install directories." FORCE)
set(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}/share")
message(STATUS "Installation directory is ${CMAKE_INSTALL_PREFIX}")

# Configure external libraries
set(EXTERNAL_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external)
add_subdirectory(${EXTERNAL_SRC_DIR})

# Configure acados
add_subdirectory(acados)

# Add examples
add_subdirectory(examples)

# Configure SWIG interface
if(NOT CMAKE_BUILD_TYPE MATCHES NoExternalLibs)
    add_subdirectory(swig)
endif()

# Generate test data and configure unit tests
if(UNIT_TESTS)
    add_subdirectory(test)
endif()

# Cpplint
include(lint)
