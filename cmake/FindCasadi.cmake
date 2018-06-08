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

# Find Casadi libraries
find_library(Casadi_LIBRARY
    NAMES casadi
    PATHS
        ENV CASADIPATH
        ${CMAKE_SOURCE_DIR}/external/*
    PATH_SUFFIXES casadi lib casadi/lib)

find_path(Casadi_INCLUDE_DIR
    NAMES casadi/casadi.hpp
    PATHS
        ENV CASADIPATH
        ${CMAKE_SOURCE_DIR}/external/*
    PATH_SUFFIXES casadi include casadi/include)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Casadi DEFAULT_MSG Casadi_LIBRARY Casadi_INCLUDE_DIR)

add_library(casadi UNKNOWN IMPORTED)
set_target_properties(casadi
    PROPERTIES
        IMPORTED_LOCATION "${Casadi_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${Casadi_INCLUDE_DIR}")

mark_as_advanced(Casadi_INCLUDE_DIR Casadi_LIBRARY)
