#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.


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
find_package_handle_standard_args(FortranLibs FOUND_VAR FORTRANLIBS_FOUND
                                              REQUIRED_VARS FORTRAN_LIBRARY)
