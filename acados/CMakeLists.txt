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

# Only test acados library for coverage
if(COVERAGE MATCHES "lcov")
    include(CodeCoverage)
    setup_target_for_coverage(acados_coverage ctest coverage)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -fprofile-arcs -ftest-coverage")
endif()

# Build list with all source files to go into the acados library
file(GLOB_RECURSE ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/*.c")
# Exclude helper files
list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/ocp_qp/condensing_helper_functions.c")
# Exclude files not compatible with MSVC
if(CMAKE_C_COMPILER_ID MATCHES "MSVC")
    list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/ocp_qp/ocp_qp_hpmpc.c")
endif()

if(NOT EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
    list(REMOVE_ITEM ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/ocp_qp/ocp_qp_ooqp.c")
endif()

# Define acados library
add_library(acados ${ACADOS_SRC})
target_include_directories(acados
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..> # for include "acados/..."
        $<BUILD_INTERFACE:${EXTERNAL_SRC_DIR}>
        $<INSTALL_INTERFACE:include>)

# HPMPC must come before BLASFEO!
if(NOT CMAKE_C_COMPILER_ID MATCHES "MSVC")
    target_link_libraries(acados PUBLIC hpmpc)
endif()
target_link_libraries(acados PUBLIC qpOASES_e blasfeo qpdunes)

# TODO: The following definition should go ASAP!
target_compile_definitions(acados PRIVATE EXT_DEP)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    target_compile_definitions(acados PRIVATE DEBUG)
endif()

find_package(OpenBLAS REQUIRED)
add_library(openblas UNKNOWN IMPORTED)
set_property(TARGET openblas PROPERTY IMPORTED_LOCATION ${OpenBLAS_LIB})

find_package(FortranLibs REQUIRED)
add_library(gfortran UNKNOWN IMPORTED)
set_property(TARGET gfortran PROPERTY IMPORTED_LOCATION ${FORTRAN_LIBRARY})

if(EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP AND NOT CMAKE_BUILD_TYPE MATCHES NoExternalLibs)
    target_compile_definitions(acados PRIVATE OOQP)
    target_link_libraries(acados
        PUBLIC
            ooqpgensparse
            ooqpsparse
            ooqpgondzio
            ooqpbase
            ma27
            openblas
            gfortran
            m)
endif()

if(COVERAGE MATCHES "lcov")
    target_link_libraries(acados PRIVATE --coverage)
endif()

target_compile_definitions(acados PRIVATE USE_ACADOS_TYPES) # Needed for qpDUNES
