##
##     This file is part of qpOASES.
##
##     qpOASES -- An Implementation of the Online Active Set Strategy.
##     Copyright (C) 2007-2017 by Hans Joachim Ferreau, Andreas Potschka,
##     Christian Kirches et al. All rights reserved.
##
##     qpOASES is free software; you can redistribute it and/or
##     modify it under the terms of the GNU Lesser General Public
##     License as published by the Free Software Foundation; either
##     version 2.1 of the License, or (at your option) any later version.
##
##     qpOASES is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##     See the GNU Lesser General Public License for more details.
##
##     You should have received a copy of the GNU Lesser General Public
##     License along with qpOASES; if not, write to the Free Software
##     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##



##
##     Filename:  CMakeLists.txt
##     Author:    Hans Joachim Ferreau (thanks to Milan Vukov)
##     Version:   3.2
##     Date:      2007-2017
##

cmake_minimum_required(VERSION 2.6)

PROJECT(qpOASES C)
SET(PACKAGE_NAME "qpOASES")
SET(PACKAGE_VERSION "3.1embeded")
SET(PACKAGE_SO_VERSION "3.1")
SET(PACKAGE_DESCRIPTION "An implementation of the online active set strategy")
SET(PACKAGE_AUTHOR "Hans Joachim Ferreau, Andreas Potschka, Christian Kirches et al.")
SET(PACKAGE_MAINTAINER "Hans Joachim Ferreau, Andreas Potschka, Christian Kirches et al.")
SET(PACKAGE_URL "https://projects.coin-or.org/qpOASES")

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

IF( NOT CMAKE_VERBOSE_MAKEFILE )
    SET( CMAKE_VERBOSE_MAKEFILE ON )
ENDIF( NOT CMAKE_VERBOSE_MAKEFILE )

IF( NOT CMAKE_BUILD_TYPE )
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
        )
ENDIF( NOT CMAKE_BUILD_TYPE )


option(QPOASES_BUILD_EXAMPLES "Build examples." ON)


############################################################
#################### compiler flags ########################
############################################################
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__NO_COPYRIGHT__")
# IF ( UNIX )
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w -O3 -finline-functions -fPIC")
    # SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_DEBUG} -O3 -finline-functions")
    ADD_DEFINITIONS(-DLINUX)
    ADD_DEFINITIONS(-D__SUPPRESSANYOUTPUT__)
    ADD_DEFINITIONS(-D__NO_STATIC__)
# ELSEIF( WINDOWS )
    # SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -nologo -EHsc -DWIN32")
# ENDIF()

# SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D__DEBUG__")

############################################################
######################## rpath #############################
############################################################
# use, i.e. don't skip the full RPATH for the build tree
# set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# # when building, don't use the install RPATH already
# # (but later on when installing)
# set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# # add the automatically determined parts of the RPATH
# # which point to directories outside the build tree to the install RPATH
# set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# # the RPATH to be used when installing, but only if it's not a system directory
# list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
# if("${isSystemDir}" STREQUAL "-1")
#    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# endif("${isSystemDir}" STREQUAL "-1")



############################################################
#################### build and install #####################
############################################################
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)

# compile qpOASES libraries
FILE(GLOB SRC src/*.c)

# library
ADD_LIBRARY(qpOASES_e STATIC ${SRC})
TARGET_LINK_LIBRARIES(qpOASES_e m)
INSTALL(TARGETS qpOASES_e
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION lib
)
SET_TARGET_PROPERTIES(
    qpOASES_e
    PROPERTIES
    SOVERSION ${PACKAGE_SO_VERSION}
    )

# headers
INSTALL(FILES include/qpOASES_e.h
  DESTINATION include)
INSTALL(DIRECTORY include/qpOASES_e
  DESTINATION include
  FILES_MATCHING PATTERN "*.h"
  PATTERN ".svn" EXCLUDE)

############################################################
######################### examples #########################
############################################################

if (QPOASES_BUILD_EXAMPLES)
    # compile qpOASES examples
    SET(EXAMPLE_NAMES
        example1
        example1b
        example3
        exampleLP
    )

    FOREACH(ELEMENT ${EXAMPLE_NAMES})
        ADD_EXECUTABLE(${ELEMENT} examples/${ELEMENT}.c)
        TARGET_LINK_LIBRARIES(${ELEMENT} qpOASES_e)
    ENDFOREACH(ELEMENT ${EXAMPLE_NAMES})
endif(QPOASES_BUILD_EXAMPLES)

##
##   end of file
##
