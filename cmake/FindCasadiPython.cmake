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
    CASADI_PYTHON_ROOT
    NAMES "casadi/casadi.py"
    HINTS ${CMAKE_SOURCE_DIR}/external/*
    PATHS ENV PYTHONPATH)

if(NOT CASADI_PYTHON_ROOT)
    message(FATAL_ERROR "Casadi not found!")
endif()
message(STATUS "Found Casadi-python: ${CASADI_PYTHON_ROOT}")

# Determine the version number
execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c "from sys import path; path.insert(0, r'${CASADI_PYTHON_ROOT}');import casadi; print(casadi.CasadiMeta_getVersion())"
    OUTPUT_VARIABLE CASADI_PYTHON_VERSION)

string(STRIP "${CASADI_PYTHON_VERSION}" CASADI_PYTHON_VERSION)
string(SUBSTRING "${CASADI_PYTHON_VERSION}" 0 1 CASADI_PYTHON_MAJOR_VERSION)
string(COMPARE EQUAL "${CASADI_PYTHON_MAJOR_VERSION}" "3" FOUND_CASADI_PYTHON_3)
if(NOT FOUND_CASADI_PYTHON_3)
    message(FATAL_ERROR "Casadi version 3 required. Found version: ${CASADI_PYTHON_VERSION}")
endif()

# Find Casadi libraries
find_library(CASADI_PYTHON_LIBRARY
    NAMES casadi
    PATHS
        "${CASADI_PYTHON_ROOT}/casadi/"
        "${CASADI_PYTHON_ROOT}/lib/"
        "${CASADI_PYTHON_ROOT}/../*")

find_path(CASADI_PYTHON_INCLUDE_DIR
    NAMES casadi/casadi.hpp
    PATHS
        "${CASADI_PYTHON_ROOT}/include"
        "${CASADI_PYTHON_ROOT}/../include")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CASADIPYTHON DEFAULT_MSG CASADI_PYTHON_LIBRARY CASADI_PYTHON_INCLUDE_DIR)

mark_as_advanced(CASADI_PYTHON_INCLUDE_DIR)
