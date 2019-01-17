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

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-error -w")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error -w")
endif()

# TODO: The following definition should go ASAP!
set(TARGET "GENERIC" CACHE STRING "Target architecture")
set(LA "HIGH_PERFORMANCE" CACHE STRING "Linear algebra optimization level")
set(BLASFEO_HEADERS_INSTALLATION_DIRECTORY "include/blasfeo/include" CACHE STRING "")
set(HPIPM_HEADERS_INSTALLATION_DIRECTORY "include/hpipm/include" CACHE STRING "")
set(HPMPC_HEADERS_INSTALLATION_DIRECTORY "include/hpmpc/include" CACHE STRING "")


if(CMAKE_SYSTEM_NAME MATCHES "dSpace")
    set(USE_C99_MATH OFF CACHE BOOL "fmax is not implemented for dSpace.")
endif()

add_subdirectory(blasfeo)
add_subdirectory(hpipm)

if(ACADOS_WITH_HPMPC)
	add_subdirectory(hpmpc)
endif()

if(ACADOS_WITH_QORE)
    set(QORE_DENSE TRUE CACHE BOOL "Build dense qore library")
    add_subdirectory(qore)
endif()

if(ACADOS_WITH_QPDUNES)
    add_subdirectory(qpdunes)
endif()

if(ACADOS_WITH_QPOASES)
    add_subdirectory(qpoases)
endif()

if(ACADOS_WITH_OSQP)
    set(PROFILING OFF CACHE BOOL "Turn off profiling in OSQP")
    set(CTRLC OFF CACHE BOOL "Turn off CTRLC in OSQP")
    add_subdirectory(osqp)
endif()

if (ACADOS_WITH_OOQP)
    include(external/ooqp)
endif()
