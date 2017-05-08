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

find_package(Matlab REQUIRED)

# Try to find Casadi root directory
find_path(
    CASADIMATLAB_ROOT_DIR
    NAMES "+casadi"
    PATHS ENV MATLABPATH
)
if (NOT CASADIMATLAB_ROOT_DIR)
    message(FATAL_ERROR "Casadi not found!")
endif ()

# Determine the version number
execute_process(
    COMMAND "${MATLAB_EXECUTABLE}" -nodesktop -nosplash -r "try, import casadi.*, disp(['casadi=',casadi.CasadiMeta.getVersion]), catch ME, disp(ME.getReport()), exit(1), end, exit(0)"
    OUTPUT_VARIABLE MATLAB_OUTPUT
)
string(FIND "${MATLAB_OUTPUT}" "casadi=" VERSION_POSITION)
string(SUBSTRING "${MATLAB_OUTPUT}" ${VERSION_POSITION} 8 CASADIMATLAB_RAW_VERSION)
string(SUBSTRING "${CASADIMATLAB_RAW_VERSION}" 7 1 CASADIMATLAB_MAJOR_VERSION)
string(COMPARE EQUAL "${CASADIMATLAB_MAJOR_VERSION}" "3" FOUND_CASADIMATLAB_3)
if (NOT FOUND_CASADIMATLAB_3)
    message(FATAL_ERROR "Casadi version 3 required. Found version: ${CASADIMATLAB_VERSION}")
endif ()

# Find Casadi libraries
set(CASADIMATLAB_LIBRARY "${CASADIMATLAB_ROOT_DIR}/lib")
find_path(
    CASADIMATLAB_INCLUDE_DIR 
    NAMES casadi/casadi.hpp
    PATHS ${CASADIMATLAB_ROOT_DIR}/include
)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CasadiMatlab DEFAULT_MSG CASADIMATLAB_LIBRARY CASADIMATLAB_INCLUDE_DIR)
