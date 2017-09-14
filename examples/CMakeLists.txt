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

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(WERROR MATCHES OFF)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -fdiagnostics-show-option")
    else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Werror -fdiagnostics-show-option")
    endif()
endif()

add_subdirectory(c)

if(SWIG_PYTHON)
    add_subdirectory(python)
endif()

if(SWIG_MATLAB)
    add_subdirectory(matlab)
endif()
