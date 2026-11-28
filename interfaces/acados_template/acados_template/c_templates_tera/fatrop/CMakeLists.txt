#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.
#

# fatrop backend for an acados OCP solver, built by AcadosFatropOcpSolver.
cmake_minimum_required(VERSION 3.16)
project(acados_fatrop C CXX)
set(CMAKE_CXX_STANDARD 17)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(ACADOS_INCLUDE_PATH "" CACHE PATH "acados include directory")
set(ACADOS_LIB_PATH "" CACHE PATH "acados lib directory")

# fatrop and acados share acados' BLASFEO
set(blasfeo_DIR ${ACADOS_LIB_PATH}/../share/cmake/blasfeo CACHE PATH "BLASFEO of acados")
find_package(fatrop REQUIRED)

add_library(acados_fatrop SHARED acados_fatrop_adapter.cpp)
target_include_directories(acados_fatrop PRIVATE
    ${ACADOS_INCLUDE_PATH}
    ${ACADOS_INCLUDE_PATH}/acados
    ${ACADOS_INCLUDE_PATH}/blasfeo/include
    ${ACADOS_INCLUDE_PATH}/hpipm/include)
target_link_directories(acados_fatrop PRIVATE ${ACADOS_LIB_PATH})
target_link_libraries(acados_fatrop PRIVATE fatrop::fatrop acados blasfeo)
set_target_properties(acados_fatrop PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    BUILD_RPATH ${ACADOS_LIB_PATH}
    # export only the C API
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # keep the symbols of a statically linked fatrop local
    target_link_options(acados_fatrop PRIVATE -Wl,--exclude-libs,ALL)
endif()
