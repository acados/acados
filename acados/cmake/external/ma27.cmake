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

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    # Needed for cross-compiling
    set(HOST_FLAG "--host=${COMPILER_PREFIX}")
endif()

ExternalProject_Add(
    ma27_project
    PREFIX "${PROJECT_BINARY_DIR}/external/ma27/"
    DOWNLOAD_COMMAND ""
    CONFIGURE_COMMAND sh "${EXTERNAL_SRC_DIR}/coinhsl/configure" "-q" "--prefix=${PROJECT_BINARY_DIR}/external/ma27/" "${HOST_FLAG}" "CC=${CMAKE_C_COMPILER}" "CFLAGS=-O2 -fPIC" "FCFLAGS=-O2 -fPIC"
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/coinhsl"
    BUILD_COMMAND make clean all
    INSTALL_COMMAND make install
    # LOG_CONFIGURE 1  # suppress output
    # LOG_BUILD 1
)

ExternalProject_Add_Step(ma27_project rename_library
    COMMAND mv ${PROJECT_BINARY_DIR}/external/ma27/lib/libcoinhsl.a ${PROJECT_BINARY_DIR}/external/ma27/lib/libma27.a
    DEPENDEES install
)

add_library(ma27 STATIC IMPORTED GLOBAL)
add_dependencies(ma27 ma27_project)
set_property(TARGET ma27 PROPERTY IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/external/ma27/lib/libma27.a")
install(FILES "${PROJECT_BINARY_DIR}/external/ma27/lib/libma27.a" DESTINATION lib)
