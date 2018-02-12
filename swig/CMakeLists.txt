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

cmake_minimum_required(VERSION 3.0)

find_package(SWIG REQUIRED)
include(UseSWIGMOD)

# Dependencies introduce some warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-field-initializers -Wno-deprecated-register") # for SWIG
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-mismatched-tags") # for CasADi
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if (CMAKE_C_COMPILER_VERSION VERSION_GREATER 5)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-maybe-uninitialized -Wno-missing-field-initializers -Wno-sign-compare")
endif()

# SWIG does not support target_include_directories()
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/external")
include_directories("${PROJECT_SOURCE_DIR}/interfaces")
# SWIG typemaps debug flag
# set(CMAKE_SWIG_FLAGS "-debug-tmused")

if(SWIG_MATLAB)
    add_subdirectory(matlab)
endif()

if(SWIG_PYTHON)
    add_subdirectory(python)
endif()
