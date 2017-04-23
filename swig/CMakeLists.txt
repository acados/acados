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

cmake_minimum_required(VERSION 3.2)

find_package(SWIG REQUIRED)

# SWIG wrapper file introduces some warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-field-initializers -Wno-deprecated-register")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-field-initializers")
endif()

include(${SWIG_USE_FILE})

# Python
find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
find_package(NumPy REQUIRED)
include_directories(${PYTHON_NUMPY_INCLUDE_DIR})

find_package(Casadi)

# SWIG does not support target_include_directories()
include_directories("${PROJECT_SOURCE_DIR}")
set_source_files_properties(acados.i PROPERTIES CPLUSPLUS ON)
# SWIG typemaps debug flag
# set(CMAKE_SWIG_FLAGS "-debug-tmused")
swig_add_module(acados python acados.i)
swig_link_libraries(acados ${PYTHON_LIBRARIES} acados)

install(TARGETS ${SWIG_MODULE_acados_REAL_NAME} DESTINATION ${LIBRARY_OUTPUT_PATH})
install(FILES ${PROJECT_BINARY_DIR}/swig/acados.py DESTINATION ${LIBRARY_OUTPUT_PATH})
install(FILES list_of_arrays.py DESTINATION ${LIBRARY_OUTPUT_PATH})

add_test(NAME ocp_nlp_py COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_BINARY_DIR}/examples ${PYTHON_EXECUTABLE} ${EXAMPLES_SRC_DIR}/python/ocp_nlp.py)
add_test(NAME ocp_qp_py COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_BINARY_DIR}/examples ${PYTHON_EXECUTABLE} ${EXAMPLES_SRC_DIR}/python/ocp_qp.py)
