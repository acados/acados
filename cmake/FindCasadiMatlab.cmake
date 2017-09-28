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
    CASADI_MATLAB_ROOT
    NAMES "+casadi"
    PATHS ENV MATLABPATH)

if(NOT CASADI_MATLAB_ROOT)
    message(FATAL_ERROR "Casadi not found!")
endif()

# Determine the version number if not in cache
if(NOT CASADI_MATLAB_MAJOR_VERSION)
    execute_process(
        COMMAND "${MATLAB_EXECUTABLE}" -nodesktop -nosplash -r "try, import casadi.*, disp(['casadi=',casadi.CasadiMeta.getVersion]), catch ME, disp(ME.getReport()), exit(1), end, exit(0)"
        OUTPUT_VARIABLE MATLAB_OUTPUT
    )
    string(FIND "${MATLAB_OUTPUT}" "casadi=" VERSION_POSITION)
    string(SUBSTRING "${MATLAB_OUTPUT}" ${VERSION_POSITION} 8 CASADI_MATLAB_RAW_VERSION)
    string(SUBSTRING "${CASADI_MATLAB_RAW_VERSION}" 7 1 CASADI_MATLAB_MAJOR_VERSION)
    set(CASADI_MATLAB_MAJOR_VERSION "${CASADI_MATLAB_MAJOR_VERSION}" CACHE STRING "")
endif()
string(COMPARE EQUAL "${CASADI_MATLAB_MAJOR_VERSION}" "3" FOUND_CASADI_MATLAB_3)
if(NOT FOUND_CASADI_MATLAB_3)
    message(FATAL_ERROR "Casadi version 3 required. Found version: ${CASADI_MATLAB_VERSION}")
endif()

# Find Casadi libraries
find_library(CASADI_MATLAB_LIBRARY
    NAMES casadi
    PATHS
        "${CASADI_MATLAB_ROOT}"
        "${CASADI_MATLAB_ROOT}/lib/"
        "${CASADI_MATLAB_ROOT}/../*")

find_path(CASADI_MATLAB_INCLUDE_DIR
    NAMES casadi/casadi.hpp
    PATHS
        "${CASADI_MATLAB_ROOT}/include"
        "${CASADI_MATLAB_ROOT}/../include")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CASADIMATLAB DEFAULT_MSG CASADI_MATLAB_LIBRARY CASADI_MATLAB_INCLUDE_DIR)

mark_as_advanced(
    CASADI_MATLAB_MAJOR_VERSION
    CASADI_MATLAB_INCLUDE_DIR
)
