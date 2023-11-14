#
# Copyright (c) The acados authors.
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

set(ACADOS_INSTALL_DIR "${PROJECT_SOURCE_DIR}" CACHE PATH  "Installation path to PROJECT_SOURCE_DIR")
if ("${CMAKE_HOST_SYSTEM_PROCESSOR}" STREQUAL "aarch64")
  set(BLASFEO_TARGET "ARMV8A_ARM_CORTEX_A57" CACHE STRING "BLASFEO Target architecture")
else()
  set(BLASFEO_TARGET "X64_AUTOMATIC" CACHE STRING "BLASFEO Target architecture")
endif()
set(HPIPM_TARGET "X64_AUTOMATIC" CACHE STRING "HPIPM Target architecture")
set(LA "HIGH_PERFORMANCE" CACHE STRING "Linear algebra optimization level")

if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(BUILD_SHARED_LIBS OFF CACHE STRING "Build shared libraries")
else()
    set(BUILD_SHARED_LIBS ON CACHE STRING "Build shared libraries")
endif()

option(ACADOS_WITH_OPENMP "OpenMP Parallelization" OFF)
option(ACADOS_SILENT "No console status output" OFF)
option(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE "Print QP inputs and outputs to file in SQP" OFF)

# Additional targets
option(ACADOS_UNIT_TESTS "Compile Unit tests" OFF)
option(ACADOS_EXAMPLES "Compile Examples" OFF)
option(ACADOS_LINT "Compile Lint" OFF)
# External libs
option(ACADOS_WITH_QPOASES "qpOASES solver" OFF)
option(ACADOS_WITH_DAQP "DAQP solver" OFF)
option(ACADOS_WITH_HPMPC "HPMPC solver" OFF)
option(ACADOS_WITH_QORE "QORE solver" OFF)
option(ACADOS_WITH_OOQP "OOQP solver" OFF)
option(ACADOS_WITH_QPDUNES "qpDUNES solver" OFF)
option(ACADOS_WITH_OSQP "OSQP solver" OFF)
# Interfaces
option(ACADOS_MATLAB "The Matlab Interface" OFF)
option(ACADOS_OCTAVE "The Octave Interface" OFF)
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

message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")

enable_testing(true)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND NOT CMAKE_SYSTEM_NAME MATCHES "Windows")
    add_definitions(-D _GLIBCXX_USE_CXX11_ABI=0) # Needed for CasADi compatibility
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_SYSTEM_NAME MATCHES "Windows")
    # MinGW: remove prefix and change suffix to match MSVC
    # (such that Matlab mex recognizes the libraries)
    set(CMAKE_SHARED_LIBRARY_PREFIX "")
    set(CMAKE_IMPORT_LIBRARY_SUFFIX ".lib")
    set(CMAKE_IMPORT_LIBRARY_PREFIX "")
    set(CMAKE_STATIC_LIBRARY_SUFFIX ".lib")
    set(CMAKE_STATIC_LIBRARY_PREFIX "")
    if(BUILD_SHARED_LIBS MATCHES ON)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wl,--export-all-symbols")
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -fdiagnostics-show-option")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -fdiagnostics-show-option")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
    if(BUILD_SHARED_LIBS MATCHES ON)
        set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS 1)
    endif()
endif()

# OPENMP
if(ACADOS_WITH_OPENMP)
    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.9.0")
        find_package(OpenMP)
    else()
        if(NOT TARGET OpenMP::OpenMP_CXX)
            set(OpenMP_CXX_FLAGS "-fopenmp")
            set(OpenMP_C_FLAGS "-fopenmp")
            find_package(Threads REQUIRED)
            add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
            set_property(TARGET OpenMP::OpenMP_CXX
                         PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
            # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
            set_property(TARGET OpenMP::OpenMP_CXX
                         PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)
            set(OpenMP_CXX_FOUND TRUE)
            set(OpenMP_C_FOUND TRUE)
        endif()
    endif()

    if(OpenMP_CXX_FOUND AND OpenMP_C_FOUND)
        message(STATUS "OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}, OpenMP_C_FLAGS: ${OpenMP_C_FLAGS}")
        # add openmp compiler flags
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

        if(ACADOS_NUM_THREADS)
            add_definitions(-DACADOS_NUM_THREADS=${ACADOS_NUM_THREADS}) # supported by CMAKE_VERSION < 3.12.0
        endif()
    else()
        message(STATUS "OpenMP NOT found.")
        set(ACADOS_WITH_OPENMP OFF)
    endif()
else()
    set(ACADOS_WITH_OPENMP OFF)
    message(STATUS "ACADOS_WITH_OPENMP: ${ACADOS_WITH_OPENMP}")
endif()

if(ACADOS_SILENT)
    message(STATUS "ACADOS_SILENT is ON")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DACADOS_SILENT")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DACADOS_SILENT")
endif()

if(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
    message(STATUS "ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE is ON")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE")
endif()

# uninstall
if(NOT TARGET uninstall)
    # Configure Uninstall
    # https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake")
endif()

# Rpath handling, see https://cmake.org/Wiki/CMake_RPATH_handling#What_is_RPATH_.3F
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

if(CMAKE_BUILD_TYPE MATCHES WithExternalLibs)
    set(ACADOS_WITH_HPMPC ON CACHE BOOL "Add HPMPC solver")
    set(ACADOS_WITH_QORE ON CACHE BOOL "Add QORE solver")
    set(ACADOS_WITH_OOQP ON CACHE BOOL "Add OOQP solver")
    set(ACADOS_WITH_QPOASES ON CACHE BOOL "Add qpOASES solver")
	set(ACADOS_WITH_DAQP ON CACHE BOOL "Add DAQP solver")
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

if(ACADOS_MATLAB OR ACADOS_OCTAVE OR ACADOS_OCTAVE_TEMPLATE OR ACADOS_PYTHON)
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

# Current config messages
message(STATUS " ")
message(STATUS "Target: BLASFEO is ${BLASFEO_TARGET}, HPIPM is ${HPIPM_TARGET}")
message(STATUS "Linear algebra: ${LA}")
message(STATUS "Matlab MEX (${ACADOS_MATLAB})")
message(STATUS "Octave MEX (${ACADOS_OCTAVE})")
message(STATUS "Octave Templates (${ACADOS_OCTAVE_TEMPLATE})")
message(STATUS "System name:version ${CMAKE_SYSTEM_NAME}:${CMAKE_SYSTEM_VERSION}")
message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")
message(STATUS "Installation directory is ${CMAKE_INSTALL_PREFIX}")
if(ACADOS_WITH_OPENMP)
    message(STATUS "OpenMP parallelization is ON")
    if(ACADOS_NUM_THREADS)
        message(STATUS "Number of threads for acados with openMP (ACADOS_NUM_THREADS) ${ACADOS_NUM_THREADS}")
    else()
        message(STATUS "Number of threads for acados with openMP (ACADOS_NUM_THREADS) not set, will use omp_get_max_threads")
    endif()
else()
    message(STATUS "OpenMP parallelization is OFF")
endif()

message(STATUS " ")

# Write libraries to be linked into a json file
if(${ACADOS_WITH_OPENMP})
    set(LINK_FLAG_OPENMP ${OpenMP_C_FLAGS})
endif()
if(${ACADOS_WITH_QPOASES})
    set(LINK_FLAG_QPOASES -lqpOASES_e)
endif()
if(${ACADOS_WITH_DAQP})
    set(LINK_FLAG_DAQP -ldaqp)
endif()
if(${ACADOS_WITH_QPDUNES})
    set(LINK_FLAG_QPDUNES -lqpdunes)
endif()
# TODO: add -lqdldl (?!)
if(${ACADOS_WITH_OSQP})
    set(LINK_FLAG_OSQP -losqp)
endif()
if(${ACADOS_WITH_HPMPC})
    set(LINK_FLAG_HPMPC -lhpmpc)
endif()
# if(${ACADOS_WITH_QORE})
#     set(LINK_FLAG_QORE -lqore)
# endif()
if(${ACADOS_WITH_OOQP})
    set(LINK_FLAG_OOQP -looqp)
endif()

file(WRITE lib/link_libs.json {\n)
file(APPEND lib/link_libs.json \t\"openmp\":\ \"${LINK_FLAG_OPENMP}\",\n)
file(APPEND lib/link_libs.json \t\"qpoases\":\ \"${LINK_FLAG_QPOASES}\",\n)
file(APPEND lib/link_libs.json \t\"daqp\":\ \"${LINK_FLAG_DAQP}\",\n)
file(APPEND lib/link_libs.json \t\"qpdunes\":\ \"${LINK_FLAG_QPDUNES}\",\n)
file(APPEND lib/link_libs.json \t\"osqp\":\ \"${LINK_FLAG_OSQP}\",\n)
file(APPEND lib/link_libs.json \t\"hpmpc\":\ \"${LINK_FLAG_HPMPC}\",\n)
# file(APPEND lib/link_libs.json \t\"qore\":\ \"${LINK_FLAG_QORE}\",\n)
file(APPEND lib/link_libs.json \t\"ooqp\":\ \"${LINK_FLAG_OOQP}\"\n) # no final comma!
file(APPEND lib/link_libs.json })
