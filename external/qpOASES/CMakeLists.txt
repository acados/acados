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
##     Author:    Hans Joachim Ferreau (thanks to Milan Vukov), Robin Verschueren
##     Version:   3.1embedded
##     Date:      2007-2017
##

cmake_minimum_required(VERSION 2.6)

project(qpOASES C)
set(PACKAGE_NAME "qpOASES")
set(PACKAGE_VERSION "3.1embeded")
set(PACKAGE_SO_VERSION "3.1")
set(PACKAGE_DESCRIPTION "An implementation of the online active set strategy")
set(PACKAGE_AUTHOR "Hans Joachim Ferreau, Andreas Potschka, Christian Kirches et al.")
set(PACKAGE_MAINTAINER "Hans Joachim Ferreau, Andreas Potschka, Christian Kirches et al.")
set(PACKAGE_URL "https://projects.coin-or.org/qpOASES")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

if(NOT CMAKE_VERBOSE_MAKEFILE)
    set(CMAKE_VERBOSE_MAKEFILE OFF)
endif(NOT CMAKE_VERBOSE_MAKEFILE)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
endif(NOT CMAKE_BUILD_TYPE)


option(QPOASES_BUILD_EXAMPLES "Build examples." ON)


############################################################
#################### compiler flags ########################
############################################################
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__NO_COPYRIGHT__")
if(CMAKE_C_COMPILER_ID MATCHES "GNU" OR CMAKE_C_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w -fPIC")
    add_definitions(-DLINUX)
    add_definitions(-D__SUPPRESSANYOUTPUT__)
    add_definitions(-D__NO_STATIC__)
elseif(CMAKE_C_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -nologo -EHsc")
    add_definitions(-DWIN32)
    add_definitions(-D__SUPPRESSANYOUTPUT__)
    add_definitions(-D__NO_STATIC__)
endif()

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D__DEBUG__")

############################################################
#################### build and install #####################
############################################################

# compile qpOASES libraries
file(GLOB SRC src/*.c)

# library
add_library(qpOASES_e ${SRC})
target_include_directories(qpOASES_e PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

if(CMAKE_C_COMPILER_ID MATCHES "GNU" OR CMAKE_C_COMPILER_ID MATCHES "Clang")
    target_link_libraries(qpOASES_e m)
endif()

install(TARGETS qpOASES_e EXPORT qpOASES_eConfig
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION lib)
install(EXPORT qpOASES_eConfig DESTINATION share)
set_target_properties(
    qpOASES_e
    PROPERTIES
    SOVERSION ${PACKAGE_SO_VERSION})
export(EXPORT qpOASES_eConfig FILE ${PROJECT_BINARY_DIR}/qpOASES_eConfig.cmake)

# headers
install(FILES include/qpOASES_e.h
  DESTINATION include)
install(DIRECTORY include/qpOASES_e
  DESTINATION include
  FILES_MATCHING PATTERN "*.h"
  PATTERN ".svn" EXCLUDE)

############################################################
######################### examples #########################
############################################################

if(QPOASES_BUILD_EXAMPLES)
    # compile qpOASES examples
    set(EXAMPLE_NAMES
        example1
        example1b
        example3
        exampleLP)
    foreach(ELEMENT ${EXAMPLE_NAMES})
        add_executable(${ELEMENT} examples/${ELEMENT}.c)
        target_link_libraries(${ELEMENT} qpOASES_e)
    endforeach()
endif()
