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

list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/sim/sim_discrete_model.c")
list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/broken_src/ocp_nlp_sm_gn.c")
list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/broken_src/ocp_nlp_sqp_bkp.c")
list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/broken_src/ocp_qp_ooqp.c")

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

# HPMPC must come before BLASFEO!
if(ACADOS_WITH_HPMPC)
    target_link_libraries(acados PUBLIC hpmpc)
endif()

if(ACADOS_WITH_QORE)
    target_link_libraries(acados PUBLIC qore)
endif()

if(ACADOS_WITH_QPOASES)
    target_link_libraries(acados PUBLIC qpOASES_e)
endif()

if(ACADOS_WITH_QPDUNES)
    target_link_libraries(acados PUBLIC qpdunes)
endif()

if(ACADOS_WITH_OSQP)
    target_link_libraries(acados PUBLIC osqp)
endif()

target_link_libraries(acados PUBLIC hpipm blasfeo)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    target_compile_definitions(acados PRIVATE DEBUG)
endif()

if(ACADOS_WITH_HPMPC)
    target_compile_definitions(acados PUBLIC ACADOS_WITH_HPMPC)
endif()

if(ACADOS_WITH_QORE)
    target_compile_definitions(acados PUBLIC ACADOS_WITH_QORE)
endif()

if(ACADOS_WITH_QPOASES)
    target_compile_definitions(acados PUBLIC ACADOS_WITH_QPOASES)
    target_compile_definitions(acados PUBLIC USE_ACADOS_TYPES)
endif()

if(ACADOS_WITH_QPDUNES)
    target_compile_definitions(acados PUBLIC ACADOS_WITH_QPDUNES)
    target_compile_definitions(acados PUBLIC USE_ACADOS_TYPES)
endif()

if(ACADOS_WITH_OSQP)
    target_compile_definitions(acados PUBLIC ACADOS_WITH_OSQP)
    # TODO(dimitris): needed or not?
    # target_compile_definitions(acados PUBLIC USE_ACADOS_TYPES)
endif()

if(ACADOS_WITH_OOQP)
    target_compile_definitions(acados PUBLIC ACADOS_WITH_OOQP)
    target_link_libraries(acados PUBLIC ooqp)
endif()

# Enable or disable timings
if(NOT CMAKE_BUILD_TYPE MATCHES "Debug")
    target_compile_definitions(acados PUBLIC MEASURE_TIMINGS)
endif()

add_subdirectory(../interfaces/acados_c acados/interfaces/acados_c)
if((SWIG_PYTHON OR SWIG_MATLAB) AND NOT CMAKE_SYSTEM_NAME MATCHES "dSpace")
    find_package(Casadi REQUIRED)
    add_subdirectory(../interfaces/acados_cpp acados/interfaces/acados_cpp)
    target_link_libraries(acados PUBLIC casadi)
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
                              INSTALL_DESTINATION cmake)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/acadosConfig.cmake
        DESTINATION cmake)

# set_target_properties(acados
#     PROPERTIES SOVERSION ${PACKAGE_SO_VERSION})
# export(TARGETS acados FILE ${PROJECT_BINARY_DIR}/acadosConfig.cmake)
