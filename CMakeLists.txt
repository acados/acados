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

cmake_minimum_required(VERSION 3.0)

if(CMAKE_MAKE_PROGRAM)
    file(TO_CMAKE_PATH ${CMAKE_MAKE_PROGRAM} CMAKE_MAKE_PROGRAM)
    find_program(CMAKE_MAKE_PROGRAM ${CMAKE_MAKE_PROGRAM})
endif()

project(acados LANGUAGES C CXX)

# if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_C_COMPILER_VERSION MATCHES "5")
#     add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)
# endif()

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

if(CMAKE_SYSTEM MATCHES "Windows")
    set(CMAKE_FIND_LIBRARY_PREFIXES "lib" CACHE INTERNAL "Library prefix")
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".dll" CACHE INTERNAL "Library suffix")
endif()

# Set up build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release
        CACHE STRING "Valid build types are: \
        Release (default), None, Debug, NoExternalLibs, RelWithDebInfo, MinSizeRel." FORCE)
endif()
message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")

enable_testing(true)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Werror -fdiagnostics-show-option")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAS} -std=c++11 -Wall -Werror -fdiagnostics-show-option")
endif()

# Set install directories
if(NOT ACADOS_INSTALL_DIR)
    set(ACADOS_INSTALL_DIR "$ENV{HOME}/local")
endif()
set(CMAKE_INSTALL_PREFIX ${ACADOS_INSTALL_DIR}
    CACHE STRING "Install path prefix, prepended onto install directories." FORCE)
set(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}/share")
message(STATUS "Installation directory is ${CMAKE_INSTALL_PREFIX}")

# Rpath handling, see https://cmake.org/Wiki/CMake_RPATH_handling#What_is_RPATH_.3F
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# set(CMAKE_INSTALL_RPATH "\$ORIGIN")
# SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP AND NOT CMAKE_BUILD_TYPE MATCHES NoExternalLibs)
    set(ACADOS_WITH_OOQP ON CACHE BOOL "Add OOQP solver")
endif()

set(CMAKE_MACOSX_RPATH TRUE)
set(EXTERNAL_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external)
# Configure acados
add_subdirectory(acados)

# Configure external libraries
add_subdirectory(${EXTERNAL_SRC_DIR})

# Add examples
if(NOT CMAKE_SYSTEM_NAME MATCHES "dSpace")
    add_subdirectory(examples)
endif()

# Configure SWIG interface
if(NOT CMAKE_BUILD_TYPE MATCHES NoExternalLibs)
    if(SWIG_MATLAB OR SWIG_PYTHON)
        add_subdirectory(swig)
    endif()
endif()

# Generate test data and configure unit tests
if(UNIT_TESTS)
    add_subdirectory(test)
endif()

# Uninstallation
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
add_custom_target(uninstall
    "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake")

# Cpplint
if (NOT (CMAKE_SYSTEM_NAME MATCHES "dSpace" OR CMAKE_GENERATOR MATCHES "MinGW Makefiles" OR CMAKE_CXX_COMPILER_ID MATCHES "MSVC"))
include(lint)
else()
message(STATUS "Skipping lint (python not found).")
endif()
