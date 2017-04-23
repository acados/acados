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

execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c "import casadi; print(casadi.CasadiMeta_getVersion())"
    OUTPUT_VARIABLE CASADI_VERSION
)
string(STRIP "${CASADI_VERSION}" CASADI_VERSION)
if (NOT CASADI_VERSION)
    message(FATAL_ERROR "Casadi not found!")
endif ()
string(SUBSTRING "${CASADI_VERSION}" 0 1 CASADI_MAJOR_VERSION)
string(COMPARE EQUAL "${CASADI_MAJOR_VERSION}" "3" FOUND_CASADI_3)
if (NOT FOUND_CASADI_3)
    message(FATAL_ERROR "Casadi version 3 required. Found version: ${CASADI_VERSION}")
else ()
    execute_process(
        COMMAND "${PYTHON_EXECUTABLE}" -c "import casadi, os, inspect; print(os.path.dirname(inspect.getfile(casadi)))"
        OUTPUT_VARIABLE CASADI_LIBRARY
    )
    string(STRIP "${CASADI_LIBRARY}" CASADI_LIBRARY)
    get_filename_component(CASADI_BASE_DIR "${CASADI_LIBRARY}/../" ABSOLUTE)
    find_path(CASADI_INCLUDE_DIR NAMES casadi/casadi.hpp
        HINTS ${CASADI_BASE_DIR}
        PATHS ${CASADI_BASE_DIR}/include
    )
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(Casadi DEFAULT_MSG CASADI_LIBRARY CASADI_INCLUDE_DIR)
endif ()

