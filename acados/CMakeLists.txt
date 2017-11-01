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

if(CMAKE_C_COMPILER_ID MATCHES "MSVC")
    set(ACADOS_WITH_HPMPC OFF)
else()
    set(ACADOS_WITH_HPMPC ON)
endif()

# Build list with all source files to go into the acados library
file(GLOB_RECURSE ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/*.c")

# Exclude files not compatible with MSVC
if(NOT ACADOS_WITH_HPMPC)
    list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/ocp_qp/ocp_qp_hpmpc.c")
endif()

if(NOT ACADOS_WITH_OOQP)
    list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/ocp_qp/ocp_qp_ooqp.c")
endif()

# Define acados library
add_library(acados ${ACADOS_SRC})

target_include_directories(acados
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> # to be able to include "acados/header.h"
        $<BUILD_INTERFACE:${EXTERNAL_SRC_DIR}>
        $<INSTALL_INTERFACE:include>)

# HPMPC must come before BLASFEO!
if(ACADOS_WITH_HPMPC)
    target_link_libraries(acados PUBLIC hpmpc)
endif()
target_link_libraries(acados PUBLIC hpipm qpOASES_e blasfeo qpdunes)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    target_compile_definitions(acados PRIVATE DEBUG)
endif()

if(ACADOS_WITH_HPMPC)
    target_compile_definitions(acados PRIVATE ACADOS_WITH_HPMPC)
endif()

if(ACADOS_WITH_OOQP)
    target_compile_definitions(acados PRIVATE OOQP)
    target_link_libraries(acados PUBLIC ooqp)
endif()

# Enable or disable timings
if(NOT ACADOS_NO_TIMINGS)
    target_compile_definitions(acados PUBLIC MEASURE_TIMINGS)
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
target_compile_definitions(acados PRIVATE USE_ACADOS_TYPES) # Needed for qpDUNES

install(TARGETS acados EXPORT acadosTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../acados
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
