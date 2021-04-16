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


set(PACKAGE_SO_VERSION "0.1")

include(CMakePackageConfigHelpers)

# Build list with all source files to go into the acados library
file(GLOB_RECURSE ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/*.c")

# Exclude files not compatible with MSVC
if(NOT ACADOS_WITH_HPMPC)
    list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/ocp_qp/ocp_qp_hpmpc.c")
endif()

if(NOT ACADOS_WITH_OOQP)
    list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/ocp_qp/ocp_qp_ooqp.c")
    list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/dense_qp/dense_qp_ooqp.c")
endif()

if(NOT ACADOS_WITH_QORE)
    list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/dense_qp/dense_qp_qore.c")
endif()

if(NOT ACADOS_WITH_QPOASES)
    list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/dense_qp/dense_qp_qpoases.c")
endif()

if(NOT ACADOS_WITH_QPDUNES)
    list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/ocp_qp/ocp_qp_qpdunes.c")
endif()

if(NOT ACADOS_WITH_OSQP)
    list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/ocp_qp/ocp_qp_osqp.c")
endif()

# Define acados library
add_library(acados ${ACADOS_SRC})

target_include_directories(acados
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> # to be able to include "acados/header.h"
        $<BUILD_INTERFACE:${EXTERNAL_SRC_DIR}/qore/QPSOLVER_DENSE/include>
        $<BUILD_INTERFACE:${EXTERNAL_SRC_DIR}/qore/QPSOLVER_DENSE/source>
        $<BUILD_INTERFACE:${EXTERNAL_SRC_DIR}/qore/KKTPACK_DENSE/include>
        $<BUILD_INTERFACE:${EXTERNAL_SRC_DIR}/qore/KKTPACK_DENSE/source>
        $<BUILD_INTERFACE:${EXTERNAL_SRC_DIR}/qore/QPCORE/include>
        $<BUILD_INTERFACE:${EXTERNAL_SRC_DIR}/qore>
        $<BUILD_INTERFACE:${EXTERNAL_SRC_DIR}/hpipm/include>
        $<BUILD_INTERFACE:${EXTERNAL_SRC_DIR}>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/interfaces>
        $<INSTALL_INTERFACE:include>)

# OPENMP
# find_package(OpenMP)
if(ACADOS_WITH_OPENMP)
    target_link_libraries(acados PUBLIC OpenMP::OpenMP_CXX)

    target_compile_definitions(acados PUBLIC ACADOS_WITH_OPENMP)
endif()

# HPMPC must come before BLASFEO!
if(ACADOS_WITH_HPMPC)
    target_link_libraries(acados PUBLIC hpmpc)
    target_compile_definitions(acados PUBLIC ACADOS_WITH_HPMPC)
endif()

if(ACADOS_WITH_QORE)
    target_link_libraries(acados PUBLIC qore)

    target_compile_definitions(acados PUBLIC ACADOS_WITH_QORE)
endif()

if(ACADOS_WITH_QPOASES)
    target_link_libraries(acados PUBLIC qpOASES_e)

    target_compile_definitions(acados PUBLIC ACADOS_WITH_QPOASES)
    target_compile_definitions(acados PUBLIC USE_ACADOS_TYPES)
endif()

if(ACADOS_WITH_QPDUNES)
    target_link_libraries(acados PUBLIC qpdunes)

    target_compile_definitions(acados PUBLIC ACADOS_WITH_QPDUNES)
    target_compile_definitions(acados PUBLIC USE_ACADOS_TYPES)
endif()

if(ACADOS_WITH_OSQP)
    target_link_libraries(acados PUBLIC osqp)

    target_compile_definitions(acados PUBLIC ACADOS_WITH_OSQP)
    # TODO(dimitris): needed or not?
    # target_compile_definitions(acados PUBLIC USE_ACADOS_TYPES)
endif()

if(ACADOS_WITH_OOQP)
    target_compile_definitions(acados PUBLIC ACADOS_WITH_OOQP)
    target_link_libraries(acados PUBLIC ooqp)
endif()

target_link_libraries(acados PUBLIC hpipm blasfeo m)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    target_compile_definitions(acados PRIVATE DEBUG)
endif()

# Enable or disable timings
if(NOT CMAKE_BUILD_TYPE MATCHES "Debug")
    target_compile_definitions(acados PUBLIC MEASURE_TIMINGS)
endif()

if(EXT_DEP MATCHES OFF)
    target_compile_definitions(acados PUBLIC BLASFEO_EXT_DEP_OFF)
endif()

# Only test acados library for coverage
if(COVERAGE MATCHES "lcov")
    include(CodeCoverage)
    setup_target_for_coverage(acados_coverage ctest coverage)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -fprofile-arcs -ftest-coverage")
    target_link_libraries(acados PRIVATE --coverage)
endif()

target_compile_definitions(acados PRIVATE ALIGNMENT=64)

install(TARGETS acados EXPORT acadosTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
    PATTERN "simplified" EXCLUDE)

install(EXPORT acadosTargets DESTINATION cmake)

configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/acadosConfig.cmake.in
                              ${CMAKE_CURRENT_BINARY_DIR}/cmake/acadosConfig.cmake
                              INSTALL_DESTINATION cmake
                              PATH_VARS ACADOS_WITH_HPMPC 
                                        ACADOS_WITH_QORE
                                        ACADOS_WITH_QPOASES
                                        ACADOS_WITH_QPDUNES
                                        ACADOS_WITH_OSQP
                                        ACADOS_WITH_OOQP)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/acadosConfig.cmake
        DESTINATION cmake)

# set_target_properties(acados
#     PROPERTIES SOVERSION ${PACKAGE_SO_VERSION})
# export(TARGETS acados FILE ${PROJECT_BINARY_DIR}/acadosConfig.cmake)
