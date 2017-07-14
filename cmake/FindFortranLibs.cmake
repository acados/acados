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

set(CMAKE_FIND_LIBRARY_PREFIXES "lib")
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

find_library(FORTRAN_LIBRARY gfortran
	PATHS
	    $ENV{GFORTRAN_HOME}/lib/gcc/mingw32/6.3.0
    HINTS
        /usr/lib/gcc/x86_64-linux-gnu/*
        /usr/local/lib/gcc/*
	    $ENV{GFORTRAN_HOME}/*
        ${CMAKE_FIND_ROOT_PATH}
    CMAKE_FIND_ROOT_PATH_BOTH)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FortranLibs FOUND_VAR FORTRANLIBS_FOUND
                                         REQUIRED_VARS FORTRAN_LIBRARY)