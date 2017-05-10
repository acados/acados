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

find_package(PythonInterp 3 REQUIRED)

# Try to find Casadi root directory
find_path(
    CASADIPYTHON_ROOT_DIR
    NAMES "casadi/casadi.py"
    PATHS ENV PYTHONPATH
)
if (NOT CASADIPYTHON_ROOT_DIR)
    message(FATAL_ERROR "Casadi not found!")
endif ()

# Determine the version number
execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c "from sys import path; path.append(r'${CASADIPYTHON_ROOT_DIR}');import casadi; print(casadi.CasadiMeta_getVersion())"
    OUTPUT_VARIABLE CASADIPYTHON_VERSION
)
string(STRIP "${CASADIPYTHON_VERSION}" CASADIPYTHON_VERSION)
string(SUBSTRING "${CASADIPYTHON_VERSION}" 0 1 CASADIPYTHON_MAJOR_VERSION)
string(COMPARE EQUAL "${CASADIPYTHON_MAJOR_VERSION}" "3" FOUND_CASADIPYTHON_3)
if (NOT FOUND_CASADIPYTHON_3)
    message(FATAL_ERROR "Casadi version 3 required. Found version: ${CASADIPYTHON_VERSION}")
endif ()

# Find Casadi libraries
set(CASADIPYTHON_LIBRARY "${CASADIPYTHON_ROOT_DIR}/casadi")
find_path(
    CASADIPYTHON_INCLUDE_DIR 
    NAMES casadi/casadi.hpp
    PATHS ${CASADIPYTHON_ROOT_DIR}/include
)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CasadiPython DEFAULT_MSG CASADIPYTHON_LIBRARY CASADIPYTHON_INCLUDE_DIR)

