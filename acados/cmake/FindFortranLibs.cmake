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


if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    list(APPEND CMAKE_FIND_LIBRARY_PREFIXES "lib")
    list(APPEND CMAKE_FIND_LIBRARY_SUFFIXES ".dll")
endif()

find_library(FORTRAN_LIBRARY NAMES libgfortran.so libgfortran.dylib gfortran
    HINTS
        /usr/lib/gcc/x86_64-linux-gnu/*
        /usr/local/lib/gcc/*
        /usr/lib/gcc/arm-linux-gnueabihf/* # for Raspbian
        ${CMAKE_FIND_ROOT_PATH}
        $ENV{PATH}
    CMAKE_FIND_ROOT_PATH_BOTH)

if(NOT FORTRAN_LIBRARY)
    find_library(FORTRAN_LIBRARY gfortran-4 gfortran-3
        HINTS
            /usr/lib/gcc/x86_64-linux-gnu/*
            /usr/local/lib/gcc/*
            ${CMAKE_FIND_ROOT_PATH}
            $ENV{PATH}
        CMAKE_FIND_ROOT_PATH_BOTH)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FORTRANLIBS FOUND_VAR FORTRANLIBS_FOUND
                                              REQUIRED_VARS FORTRAN_LIBRARY)
