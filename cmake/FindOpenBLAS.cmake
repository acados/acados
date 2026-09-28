#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.


SET(Open_BLAS_INCLUDE_SEARCH_PATHS
  "/usr/include"
  "/usr/include/openblas"
  "/usr/include/openblas-base"
  "/usr/local/include"
  "/usr/local/include/openblas"
  "/usr/local/include/openblas-base"
  "/usr/local/Cellar/openblas/*/include"
  "/opt/OpenBLAS/include"
  "/opt/openblas/include"
  "$ENV{OpenBLAS_HOME}"
  "$ENV{OpenBLAS_HOME}/include"
  "$ENV{OpenBLAS_HOME}/include/OpenBLAS"
)

SET(Open_BLAS_LIB_SEARCH_PATHS
    "/lib/"
    "/lib/openblas-base"
    "/lib64/"
    "/usr/lib"
    "/usr/lib/openblas-base"
    "/usr/lib64"
    "/usr/local/lib"
    "/usr/local/lib64"
    "/usr/local/Cellar/openblas/*/lib"
    "/opt/OpenBLAS/lib"
    "/opt/openblas/lib"
    "$ENV{OpenBLAS}cd"
    "$ENV{OpenBLAS}/lib"
    "$ENV{OpenBLAS_HOME}"
    "$ENV{OpenBLAS_HOME}/lib"
    "$ENV{OpenBLAS_HOME}/bin"
    "$ENV{PATH}"
)

FIND_PATH(OpenBLAS_INCLUDE_DIR NAMES cblas.h PATHS  "$ENV{OpenBLAS_HOME}/include/OpenBLAS")
FIND_LIBRARY(OpenBLAS_LIB NAMES openblas PATHS ${Open_BLAS_LIB_SEARCH_PATHS})

SET(OpenBLAS_FOUND ON)

IF(OpenBLAS_LIB)
    get_filename_component(OpenBLAS_HOME ${OpenBLAS_LIB} DIRECTORY)
    get_filename_component(OpenBLAS_HOME ${OpenBLAS_HOME} DIRECTORY)
    FIND_PATH(OpenBLAS_INCLUDE_DIR NAMES cblas.h
        PATHS
            ${OpenBLAS_HOME}/include/OpenBLAS
            ${OpenBLAS_HOME}/include/)
ENDIF()

#    Check include files
IF(NOT OpenBLAS_INCLUDE_DIR)
    SET(OpenBLAS_FOUND OFF)
    MESSAGE(STATUS "Could not find OpenBLAS include. Turning OpenBLAS_FOUND off")
ENDIF()

#    Check libraries
IF(NOT OpenBLAS_LIB)
    SET(OpenBLAS_FOUND OFF)
    MESSAGE(STATUS "Could not find OpenBLAS lib. Turning OpenBLAS_FOUND off")
ENDIF()

IF (OpenBLAS_FOUND)
  IF (NOT OpenBLAS_FIND_QUIETLY)
    # MESSAGE(STATUS "Found OpenBLAS libraries: ${OpenBLAS_LIB}")
    # MESSAGE(STATUS "Found OpenBLAS include: ${OpenBLAS_INCLUDE_DIR}")
  ENDIF (NOT OpenBLAS_FIND_QUIETLY)
ELSE (OpenBLAS_FOUND)
  IF (OpenBLAS_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Could not find OpenBLAS")
  ENDIF (OpenBLAS_FIND_REQUIRED)
ENDIF (OpenBLAS_FOUND)

MARK_AS_ADVANCED(
    OpenBLAS_INCLUDE_DIR
    OpenBLAS_LIB
    OpenBLAS
)
