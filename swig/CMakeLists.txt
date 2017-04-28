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
include(UseSWIGMOD)

# SWIG wrapper file introduces some warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-field-initializers -Wno-deprecated-register")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-field-initializers")
endif()

# SWIG does not support target_include_directories()
include_directories("${PROJECT_SOURCE_DIR}")
set_source_files_properties(acados.i PROPERTIES CPLUSPLUS ON)
# SWIG typemaps debug flag
# set(CMAKE_SWIG_FLAGS "-debug-tmused")

set(MATLAB_SWIG_INTERFACE OFF) # Until we have MATLAB on Travis
set(PYTHON_SWIG_INTERFACE ON)
if(MATLAB_SWIG_INTERFACE)
    find_package(Matlab REQUIRED)

    swig_module_initialize(acadosMEX matlab)
    swig_add_source_to_module(acadosMEX FALSE MEX_FILE acados.i)

    if (EXISTS ${MATLAB_INCLUDE_DIR}/../version/cpp_mexapi_version.cpp)
        add_library(acadosMATLAB MODULE ${MEX_FILE} ${MATLAB_INCLUDE_DIR}/../version/cpp_mexapi_version.cpp)
    else()
        add_library(acadosMATLAB MODULE ${MEX_FILE})
    endif()
    target_include_directories(acadosMATLAB PUBLIC ${MATLAB_INCLUDE_DIR})
    target_link_libraries(acadosMATLAB ${MATLAB_LIBRARIES} acados)
    set_target_properties(acadosMATLAB
        PROPERTIES
            OUTPUT_NAME "acadosMEX"
            PREFIX ""
            SUFFIX .${MATLAB_MEX_EXT}
    )
endif()
if(PYTHON_SWIG_INTERFACE) # Python
    find_package(PythonInterp 3 REQUIRED)
    find_package(PythonLibs 3 REQUIRED)
    find_package(NumPy REQUIRED)
    find_package(Casadi)

    include_directories(${PYTHON_INCLUDE_DIRS})
    include_directories(${PYTHON_NUMPY_INCLUDE_DIR})
    swig_add_module(acados python acados.i)
    swig_link_libraries(acados ${PYTHON_LIBRARIES} acados)

    install(TARGETS ${SWIG_MODULE_acados_REAL_NAME} DESTINATION ${LIBRARY_OUTPUT_PATH})
    install(FILES ${PROJECT_BINARY_DIR}/swig/acados.py DESTINATION ${LIBRARY_OUTPUT_PATH})
    # add_test(NAME ocp_nlp_py COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_BINARY_DIR}/examples ${PYTHON_EXECUTABLE} ${EXAMPLES_SRC_DIR}/python/ocp_nlp.py)
    add_test(NAME ocp_qp_py COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_BINARY_DIR}/examples ${PYTHON_EXECUTABLE} ${EXAMPLES_SRC_DIR}/python/ocp_qp.py)
endif()
