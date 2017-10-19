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

# TODO: The following definition should go ASAP!
if(CMAKE_SYSTEM_NAME MATCHES "dSpace")
    set(EXT_DEP OFF CACHE BOOL "Compile external dependencies in BLASFEO")
endif()
set(TARGET "GENERIC" CACHE STRING "Target architecture")
set(LA "HIGH_PERFORMANCE" CACHE STRING "Linear algebra optimization level")
set(BLASFEO_HEADERS_INSTALLATION_DIRECTORY "include/blasfeo/include" CACHE STRING "")
set(HPIPM_HEADERS_INSTALLATION_DIRECTORY "include/hpipm/include" CACHE STRING "")
set(HPMPC_HEADERS_INSTALLATION_DIRECTORY "include/hpmpc/include" CACHE STRING "")

add_subdirectory(blasfeo)
add_subdirectory(hpipm)
if((NOT CMAKE_C_COMPILER_ID MATCHES "MSVC") AND (NOT CMAKE_SYSTEM_NAME MATCHES "dSpace"))
	add_subdirectory(hpmpc)
endif()
add_subdirectory(qpdunes)
add_subdirectory(qpOASES)
if (ACADOS_WITH_OOQP)
    include(external/ooqp)
endif()
