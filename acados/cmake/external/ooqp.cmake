#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#


include(ExternalProject)

find_package(OpenBLAS REQUIRED)
add_library(openblas UNKNOWN IMPORTED)
set_property(TARGET openblas PROPERTY IMPORTED_LOCATION ${OpenBLAS_LIB})

find_package(FortranLibs REQUIRED)
add_library(gfortran UNKNOWN IMPORTED)
set_property(TARGET gfortran PROPERTY IMPORTED_LOCATION ${FORTRAN_LIBRARY})

include(external/ma27)

set(OOQP_LDFLAGS "")
set(HOST_FLAG "")
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(OOQP_LDFLAGS "-lc++")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    # Needed for cross-compiling
    set(HOST_FLAG "--host=${COMPILER_PREFIX}")
endif()

ExternalProject_Add(
    ooqp_project

    DEPENDS ma27
    PREFIX "${PROJECT_BINARY_DIR}/external/OOQP"
    DOWNLOAD_COMMAND ""
    CONFIGURE_COMMAND ./configure "-q" "--prefix=${PROJECT_BINARY_DIR}/external/OOQP/" "${HOST_FLAG}"
                                  "CXX=${CMAKE_CXX_COMPILER}" "CXXFLAGS=-O2 -fPIC" "CC=${CMAKE_C_COMPILER}"
                                  "CFLAGS=-O2 -fPIC" "FFLAGS=-O2 -fPIC" "LDFLAGS=${OOQP_LDFLAGS}"
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/OOQP"
    BUILD_COMMAND make clean all
    INSTALL_COMMAND make install
    # LOG_CONFIGURE 1  # suppress output
    # LOG_BUILD 1
)

ExternalProject_Get_Property(ooqp_project BINARY_DIR)

ExternalProject_Add_Step(ooqp_project copy_ooqp
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${EXTERNAL_SRC_DIR}/OOQP/" "${BINARY_DIR}"
    DEPENDERS configure
)

ExternalProject_Add_Step(ooqp_project create_lib_folder
    COMMAND ${CMAKE_COMMAND} -E chdir ${BINARY_DIR} mkdir -p lib
    DEPENDERS configure
)

ExternalProject_Add_Step(ooqp_project copy_ma27
    COMMAND cp ${PROJECT_BINARY_DIR}/external/ma27/lib/libma27.a ${BINARY_DIR}/libma27.a
    DEPENDERS configure
)

add_library(ooqp INTERFACE)
target_link_libraries(ooqp INTERFACE
    ooqpgensparse
    ooqpsparse
    ooqpgendense
    ooqpdense
    ooqpgondzio
    ooqpbase
    ma27
    openblas
    lapack
    gfortran
    m)

set_property(TARGET ooqp
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/external/OOQP/include>
        $<INSTALL_INTERFACE:include>)

install(TARGETS ooqp EXPORT ooqpConfig
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin)

install(EXPORT ooqpConfig DESTINATION cmake)

install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/../FindFortranLibs.cmake
        ${CMAKE_CURRENT_LIST_DIR}/../FindOpenBLAS.cmake
    DESTINATION cmake)

install(DIRECTORY ${BINARY_DIR}/include/
    DESTINATION include/ooqp
    FILES_MATCHING PATTERN "*.h")

add_library(ooqpgensparse STATIC IMPORTED GLOBAL)
add_dependencies(ooqpgensparse ooqp_project)
set_property(TARGET ooqpgensparse PROPERTY IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/external/OOQP/lib/libooqpgensparse.a")
set_property(TARGET ooqpgensparse PROPERTY IMPORTED_LINK_INTERFACE_LANGUAGES CXX)
install(FILES "${PROJECT_BINARY_DIR}/external/OOQP/lib/libooqpgensparse.a" DESTINATION lib)

add_library(ooqpsparse STATIC IMPORTED GLOBAL)
add_dependencies(ooqpsparse ooqp_project)
set_property(TARGET ooqpsparse PROPERTY IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/external/OOQP/lib/libooqpsparse.a")
set_property(TARGET ooqpsparse PROPERTY IMPORTED_LINK_INTERFACE_LANGUAGES CXX)
install(FILES "${PROJECT_BINARY_DIR}/external/OOQP/lib/libooqpsparse.a" DESTINATION lib)

add_library(ooqpgendense STATIC IMPORTED GLOBAL)
add_dependencies(ooqpgendense ooqp_project)
set_property(TARGET ooqpgendense PROPERTY IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/external/OOQP/lib/libooqpgendense.a")
set_property(TARGET ooqpgendense PROPERTY IMPORTED_LINK_INTERFACE_LANGUAGES CXX)
install(FILES "${PROJECT_BINARY_DIR}/external/OOQP/lib/libooqpgendense.a" DESTINATION lib)

add_library(ooqpdense STATIC IMPORTED GLOBAL)
add_dependencies(ooqpdense ooqp_project)
set_property(TARGET ooqpdense PROPERTY IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/external/OOQP/lib/libooqpdense.a")
set_property(TARGET ooqpdense PROPERTY IMPORTED_LINK_INTERFACE_LANGUAGES CXX)
install(FILES "${PROJECT_BINARY_DIR}/external/OOQP/lib/libooqpdense.a" DESTINATION lib)

add_library(ooqpgondzio STATIC IMPORTED GLOBAL)
add_dependencies(ooqpgondzio ooqp_project)
set_property(TARGET ooqpgondzio PROPERTY IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/external/OOQP/lib/libooqpgondzio.a")
set_property(TARGET ooqpgondzio PROPERTY IMPORTED_LINK_INTERFACE_LANGUAGES CXX)
install(FILES "${PROJECT_BINARY_DIR}/external/OOQP/lib/libooqpgondzio.a" DESTINATION lib)

add_library(ooqpbase STATIC IMPORTED GLOBAL)
add_dependencies(ooqpbase ma27 ooqp_project)
set_property(TARGET ooqpbase PROPERTY IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/external/OOQP/lib/libooqpbase.a")
set_property(TARGET ooqpbase PROPERTY IMPORTED_LINK_INTERFACE_LANGUAGES CXX)
install(FILES "${PROJECT_BINARY_DIR}/external/OOQP/lib/libooqpbase.a" DESTINATION lib)
