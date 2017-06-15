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

# Dependencies introduce some warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-field-initializers -Wno-deprecated-register") # for SWIG
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-mismatched-tags") # for CasADi
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation -Wno-missing-field-initializers")
endif()

# SWIG does not support target_include_directories()
include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/external/casadi-octave-v3.1.1/include")
set_source_files_properties(acados.i PROPERTIES CPLUSPLUS ON)
# SWIG typemaps debug flag
# set(CMAKE_SWIG_FLAGS "-debug-tmused")

find_package(Matlab)
if(SWIG_MATLAB AND MATLAB_FOUND)
    find_package(Matlab REQUIRED)
    find_package(CasadiMatlab)

    if(CASADIMATLAB_FOUND)
        add_library(casadi SHARED IMPORTED)
        set_property(TARGET casadi PROPERTY IMPORTED_LOCATION "${CASADIMATLAB_LIBRARY}/libcasadi.dylib")
        install(FILES "${CASADIMATLAB_LIBRARY}/libcasadi.dylib" DESTINATION lib) 
    endif()

    swig_module_initialize(acadosMEX matlab)
    swig_add_source_to_module(acadosMEX FALSE MEX_FILE acados.i)

    if(EXISTS ${MATLAB_INCLUDE_DIR}/../version/cpp_mexapi_version.cpp)
        add_library(acadosMATLAB MODULE ${MEX_FILE} ${MATLAB_INCLUDE_DIR}/../version/cpp_mexapi_version.cpp)
    else()
        add_library(acadosMATLAB MODULE ${MEX_FILE})
    endif()
    target_include_directories(acadosMATLAB PUBLIC ${MATLAB_INCLUDE_DIR})
    target_link_libraries(acadosMATLAB ${MATLAB_LIBRARIES} acados casadi)
    set_target_properties(acadosMATLAB
        PROPERTIES
            OUTPUT_NAME "acadosMEX"
            PREFIX ""
            SUFFIX .${MATLAB_MEX_EXT}
    )
    install(TARGETS acadosMATLAB DESTINATION lib)
    install(DIRECTORY ${PROJECT_BINARY_DIR}/swig/+acados DESTINATION lib)
    install(
        FILES
            ${PROJECT_BINARY_DIR}/swig/SwigGet.m
            ${PROJECT_BINARY_DIR}/swig/SwigMem.m
            ${PROJECT_BINARY_DIR}/swig/SwigRef.m
        DESTINATION lib
    )
endif()

find_package(PythonInterp 3)
if(SWIG_PYTHON AND PYTHONINTERP_FOUND)
    find_package(PythonInterp 3 REQUIRED)
    find_package(PythonLibs 3 REQUIRED)
    find_package(NumPy REQUIRED)
    find_package(CasadiPython)

    # do not install casadi twice
    if(NOT (SWIG_MATLAB AND CASADIMATLAB_FOUND) AND CASADIPYTHON_FOUND)
        add_library(casadi SHARED IMPORTED)
        set_property(TARGET casadi PROPERTY IMPORTED_LOCATION "${CASADIPYTHON_LIBRARY}/libcasadi.dylib")
        install(FILES "${CASADIPYTHON_LIBRARY}/libcasadi.dylib" DESTINATION lib) 
    endif()

    include_directories(${PYTHON_INCLUDE_DIRS})
    include_directories(${PYTHON_NUMPY_INCLUDE_DIR})
    swig_add_module(acados python acados.i)
    swig_link_libraries(acados ${PYTHON_LIBRARIES} acados casadi)

    install(TARGETS ${SWIG_MODULE_acados_REAL_NAME} DESTINATION lib)
    install(FILES ${PROJECT_BINARY_DIR}/swig/acados.py DESTINATION lib)
endif()
