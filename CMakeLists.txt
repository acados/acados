#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#


cmake_minimum_required(VERSION 3.1)

if(CMAKE_MAKE_PROGRAM)
    file(TO_CMAKE_PATH ${CMAKE_MAKE_PROGRAM} CMAKE_MAKE_PROGRAM)
    find_program(CMAKE_MAKE_PROGRAM ${CMAKE_MAKE_PROGRAM})
endif()

project(acados LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

### Defaults

set(ACADOS_INSTALL_DIR "${PROJECT_SOURCE_DIR}" CACHE STRING  "Installation path to PROJECT_SOURCE_DIR")
set(BLASFEO_TARGET "X64_AUTOMATIC" CACHE STRING "BLASFEO Target architecture")
set(HPIPM_TARGET "AVX" CACHE STRING "HPIPM Target architecture")
set(LA "HIGH_PERFORMANCE" CACHE STRING "Linear algebra optimization level")

if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(BUILD_SHARED_LIBS OFF CACHE STRING "Build shared libraries")
else()
    set(BUILD_SHARED_LIBS ON CACHE STRING "Build shared libraries")
endif()

# Additional targets
option(ACADOS_UNIT_TESTS "Compile Unit tests" OFF)
option(ACADOS_EXAMPLES "Compile Examples" OFF)
option(ACADOS_LINT "Compile Lint" OFF)
# Extarnal libs
option(ACADOS_WITH_QPOASES  "qpOASES solver" OFF)
option(ACADOS_WITH_HPMPC "HPMPC solver" OFF)
option(ACADOS_WITH_QORE "QORE solver" OFF)
option(ACADOS_WITH_OOQP "OOQP solver" OFF)
option(ACADOS_WITH_QPDUNES "qpDUNES solver" OFF)
option(ACADOS_WITH_OSQP "OSQP solver" OFF)
# Interfaces
option(ACADOS_MATLAB "The Matlab Interface" OFF)
option(ACADOS_OCTAVE "The Octave Interface" OFF)
option(ACADOS_OCTAVE_TEMPLATE "Template MEX Interface" OFF)
option(ACADOS_PYTHON "Python Interface" OFF)

# Set custom path
set(CMAKE_INSTALL_PREFIX ${ACADOS_INSTALL_DIR})
set(EXTERNAL_SRC_DIR ${PROJECT_SOURCE_DIR}/external)
set(BLASFEO_SRC_DIR ${PROJECT_SOURCE_DIR}/external/blasfeo)

# Set up build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release
        CACHE STRING "Valid build types are: \
        Release (default), None, Debug, WithExternalLibs, RelWithDebInfo, MinSizeRel." FORCE)
endif()
### Defaults

message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")

enable_testing(true)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND NOT CMAKE_SYSTEM_NAME MATCHES "Windows")
    add_definitions(-D _GLIBCXX_USE_CXX11_ABI=0) # Needed for CasADi compatibility
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_SYSTEM_NAME MATCHES "Windows")
    # MinGW, change to .lib such that mex recognizes it
    set(CMAKE_STATIC_LIBRARY_SUFFIX ".lib")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -fdiagnostics-show-option")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fdiagnostics-show-option")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
endif()

# Rpath handling, see https://cmake.org/Wiki/CMake_RPATH_handling#What_is_RPATH_.3F
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

if(CMAKE_BUILD_TYPE MATCHES WithExternalLibs)
    set(ACADOS_WITH_HPMPC ON CACHE BOOL "Add HPMPC solver")
    set(ACADOS_WITH_QORE ON CACHE BOOL "Add QORE solver")
    set(ACADOS_WITH_OOQP ON CACHE BOOL "Add OOQP solver")
    set(ACADOS_WITH_QPOASES ON CACHE BOOL "Add qpOASES solver")
    set(ACADOS_WITH_QPDUNES ON CACHE BOOL "Add qpDUNES solver")
    set(ACADOS_WITH_OSQP ON CACHE BOOL "Add OSQP solver")
endif()

## External lib checks (

if(ACADOS_WITH_HPMPC MATCHES ON AND CMAKE_C_COMPILER_ID MATCHES "MSVC")
    set(ACADOS_WITH_HPMPC OFF CACHE BOOL "HPMPC solver" FORCE)
    message(WARNING "HPMPC has been disabled, not compatible with ${CMAKE_C_COMPILER_ID}")
endif()

if(ACADOS_WITH_OOQP MATCHES ON)

    if(CMAKE_SYSTEM_NAME MATCHES "dSpace")
        set(ACADOS_WITH_OOQP OFF CACHE BOOL "Add OOQP solver" FORCE)
        message(WARNING "OOQP is not compatible with dSpace, OOQP is disabled")

    elseif(NOT EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
        set(ACADOS_WITH_OOQP OFF CACHE BOOL "Add OOQP solver" FORCE)
        message(WARNING "OOQP folder not found, OOQP has been disabled")

    endif()
endif()

if(ACADOS_WITH_QORE MATCHES ON)

    if(CMAKE_SYSTEM_NAME MATCHES "dSpace" OR CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        set(ACADOS_WITH_QORE OFF CACHE BOOL "Add QORE solver" FORCE)
        message(WARNING "QORE is not compatible with MSVC or dSpace, QORE is has been disabled")

    elseif(NOT EXISTS ${PROJECT_SOURCE_DIR}/external/qore)
        set(ACADOS_WITH_QORE OFF CACHE BOOL "Add QORE solver" FORCE)
        message(WARNING "QORE folder not found, QORE has been disabled")

    endif()
endif()

if(ACADOS_LINT)
    if(CMAKE_SYSTEM_NAME MATCHES "Windows" OR CMAKE_SYSTEM_NAME MATCHES "dSpace")
        set(ACADOS_LINT OFF CACHE BOOL "Lint" FORCE)
        message(WARNING "Lint has been disabled, not compatible with ${CMAKE_SYSTEM_NAME}")
    endif()
endif()

if(ACADOS_MATLAB OR ACADOS_OCTAVE OR ACADOS_PYTHON) 
    add_subdirectory(${PROJECT_SOURCE_DIR}/interfaces)
endif()

## ) External lib checks

# Configure acados
add_subdirectory(acados)

# Configure interfaces
add_subdirectory(${PROJECT_SOURCE_DIR}/interfaces/acados_c)

# Configure external libraries
add_subdirectory(${EXTERNAL_SRC_DIR})

# Configure examples
if(ACADOS_EXAMPLES)
    add_subdirectory(examples)
endif()

# Configure tests
if(ACADOS_UNIT_TESTS)
    add_subdirectory(test)
endif()

# Configure lint
if(ACADOS_LINT)
    include(Lint)
else()
    add_custom_target(lint echo "Lint is disabled")
endif()

# Configure Uninstall
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
add_custom_target(uninstall
    # NOTE(dimitris): uninstall name clashes with submodules
    "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake")

# Current config messages
message(STATUS "Target: BLASFEO is ${BLASFEO_TARGET}, HPIPM is ${HPIPM_TARGET}")
message(STATUS "Linear algebra: ${LA}")
message(STATUS "Matlab MEX (${ACADOS_MATLAB})")
message(STATUS "Octave MEX (${ACADOS_OCTAVE})")
message(STATUS "System name:version ${CMAKE_SYSTEM_NAME}:${CMAKE_SYSTEM_VERSION}")
message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")
message(STATUS "Installation directory is ${CMAKE_INSTALL_PREFIX}")
