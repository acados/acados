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

find_library(FORTRAN_LIBRARY
    NAMES gfortran
    PATHS
        "/usr/local/lib/gcc/*"
        "/usr/lib/gcc/x86_64-linux-gnu/*"
        "/usr/lib/gcc/arm-linux-gnueabihf/*"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FortranLibs DEFAULT_MSG FORTRAN_LIBRARY)
