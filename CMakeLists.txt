#
#	 This file is part of acados.
#
#	 acados is free software; you can redistribute it and/or
#	 modify it under the terms of the GNU Lesser General Public
#	 License as published by the Free Software Foundation; either
#	 version 3 of the License, or (at your option) any later version.
#
#	 acados is distributed in the hope that it will be useful,
#	 but WITHOUT ANY WARRANTY; without even the implied warranty of
#	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#	 Lesser General Public License for more details.
#
#	 You should have received a copy of the GNU Lesser General Public
#	 License along with acados; if not, write to the Free Software Foundation,
#	 Inc., 51 Franklin Street, Fifth Floor, Boston, MA	02110-1301	USA
#
#

cmake_minimum_required(VERSION 3.1)

if(CMAKE_MAKE_PROGRAM)
	file(TO_CMAKE_PATH ${CMAKE_MAKE_PROGRAM} CMAKE_MAKE_PROGRAM)
	find_program(CMAKE_MAKE_PROGRAM ${CMAKE_MAKE_PROGRAM})
endif()

project(acados LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

### Defaults

set(ACADOS_INSTALL_DIR "/usr/local" CACHE STRING  "Installation path to /usr/local")
set(BLASFEO_TARGET "GENERIC" CACHE STRING "BLASFEO Target architecture")
set(HPIPM_TARGET "GENERIC" CACHE STRING "HPIPM Target architecture")
set(LA "HIGH_PERFORMANCE" CACHE STRING "Linear algebra optimization level")

# Additional targets
option(ACADOS_UNIT_TESTS "Compile Unit tests" OFF)
option(ACADOS_EXAMPLES "Compile Examples" OFF)
option(ACADOS_LINT "Compile Lint" OFF)
# Extarnal libs
option(ACADOS_WITH_QPOASES  "qpOASES solver" ON)
option(ACADOS_WITH_HPMPC "HPMPC solver" OFF)
option(ACADOS_WITH_QORE "QORE solver" OFF)
option(ACADOS_WITH_OOQP "OOQP solver" OFF)
option(ACADOS_WITH_QPDUNES "qpDUNES solver" OFF)
option(ACADOS_WITH_OSQP "OSQP solver" OFF)
# Interfaces
option(SWIG_MATLAB "Swig Matlab Interface" OFF)
option(SWIG_PYTHON "Swig Python Interface" OFF)
option(TEMPLATE_MATLAB "Template Matlab Interface" OFF)
option(TEMPLATE_PYTHON "Template Python Interface" OFF)
option(DEV_MATLAB "Developer Matlab Interface" OFF)

# Set custom path
set(CMAKE_INSTALL_PREFIX ${ACADOS_INSTALL_DIR})
set(EXTERNAL_SRC_DIR ${PROJECT_SOURCE_DIR}/external)
set(BLASFEO_SRC_DIR ${PROJECT_SOURCE_DIR}/external/blasfeo)

# Set up build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release
		CACHE STRING "Valid build types are: \
		Release (default), None, Debug, WithExternalLibs, RelWithDebInfo, MinSizeRel." FORCE)
endif()
### Defaults

message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")

enable_testing(true)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND NOT CMAKE_SYSTEM_NAME MATCHES "Windows")
	add_definitions(-D _GLIBCXX_USE_CXX11_ABI=0) # Needed for CasADi compatibility
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_SYSTEM_NAME MATCHES "Windows")
	# MinGW, change to .lib such that mex recognizes it
	set(CMAKE_STATIC_LIBRARY_SUFFIX ".lib")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -fdiagnostics-show-option")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fdiagnostics-show-option")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
	set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
endif()

# Rpath handling, see https://cmake.org/Wiki/CMake_RPATH_handling#What_is_RPATH_.3F
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "$ENV{CASADIPATH}/casadi")

if(CMAKE_BUILD_TYPE MATCHES WithExternalLibs)
	set(ACADOS_WITH_HPMPC ON CACHE BOOL "Add HPMPC solver")
	set(ACADOS_WITH_QORE ON CACHE BOOL "Add QORE solver")
	set(ACADOS_WITH_OOQP ON CACHE BOOL "Add OOQP solver")
	set(ACADOS_WITH_QPOASES ON CACHE BOOL "Add qpOASES solver")
	set(ACADOS_WITH_QPDUNES ON CACHE BOOL "Add qpDUNES solver")
	set(ACADOS_WITH_OSQP ON CACHE BOOL "Add OSQP solver")
endif()

## External lib checks (

if(ACADOS_WITH_HPMPC MATCHES ON AND CMAKE_C_COMPILER_ID MATCHES "MSVC")
	set(ACADOS_WITH_HPMPC OFF CACHE BOOL "HPMPC solver" FORCE)
	message(WARNING "HPMPC is not compatible with MSVC")
endif()

if(ACADOS_WITH_OOQP MATCHES ON)

	if(CMAKE_SYSTEM_NAME MATCHES "dSpace")
		set(ACADOS_WITH_OOQP OFF CACHE BOOL "Add OOQP solver" FORCE)
		message(WARNING "OOQP is not compatible with dSpace, OOQP is disabled")

	elseif(NOT EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
		set(ACADOS_WITH_OOQP OFF CACHE BOOL "Add OOQP solver" FORCE)
		message(WARNING "OOQP folder not found, OOQP has been disabled")

	endif()
endif()

if(ACADOS_WITH_QORE MATCHES ON)

	if(CMAKE_SYSTEM_NAME MATCHES "dSpace" OR CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
		set(ACADOS_WITH_QORE OFF CACHE BOOL "Add QORE solver" FORCE)
		message(WARNING "QORE is not compatible with MSVC or dSpace, QORE is has been disabled")

	elseif(NOT EXISTS ${PROJECT_SOURCE_DIR}/external/qore)
		set(ACADOS_WITH_QORE OFF CACHE BOOL "Add QORE solver" FORCE)
		message(WARNING "QORE folder not found, QORE has been disabled")

	endif()
endif()

if(SWIG_MATLAB OR SWIG_PYTHON)
	if(CMAKE_SYSTEM_NAME MATCHES "dSpace")
		set(SWIG_MATLAB OFF CACHE BOOL "Swig Matlab Interface" FORCE)
		set(SWIG_PYTHON OFF CACHE BOOL "Swig Python Interface" FORCE)
		message(WARNING "Swig has been disabled, not compatible with dSpace system")
	endif()
	if(CMAKE_BUILD_TYPE MATCHES NoExternalLibs)
		set(SWIG_MATLAB OFF CACHE BOOL "Swig Matlab Interface" FORCE)
		set(SWIG_PYTHON OFF CACHE BOOL "Swig Python Interface" FORCE)
		message(WARNING "Swig has been disabled, not compatible with NoExternalLibs target")
	endif()
endif()

if(ACADOS_LINT)
	if(CMAKE_SYSTEM_NAME MATCHES "Windows" OR CMAKE_SYSTEM_NAME MATCHES "dSpace")
		set(ACADOS_LINT OFF CACHE BOOL "Lint" FORCE)
		message(WARNING "Lint has been disabled, not compatible with ${CMAKE_SYSTEM_NAME}")
	endif()
endif()

### ) External lib checks

# Configure acados
add_subdirectory(acados)

# Configure interfaces
add_subdirectory(${PROJECT_SOURCE_DIR}/interfaces/acados_c)

if((SWIG_PYTHON OR SWIG_MATLAB))
    add_subdirectory(${PROJECT_SOURCE_DIR}/interfaces/acados_cpp)
    add_subdirectory(${PROJECT_SOURCE_DIR}/swig)
endif()

# Configure external libraries
add_subdirectory(${EXTERNAL_SRC_DIR})

# Configure examples
if(ACADOS_EXAMPLES)
	add_subdirectory(examples)
endif()

# Configure tests
if(ACADOS_UNIT_TESTS)
	add_subdirectory(test)
endif()

# Configure lint
if(ACADOS_LINT)
	include(Lint)
else()
	add_custom_target(lint echo "Lint is disabled")
endif()

# Configure Uninstall
configure_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
	IMMEDIATE @ONLY)
add_custom_target(uninstall
	# NOTE(dimitris): uninstall name clashes with submodules
	"${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake")

# Current config messages
message(STATUS "Target: BLASFEO is ${BLASFEO_TARGET}, HPIPM is ${HPIPM_TARGET}")
message(STATUS "Linear algebra: ${LA}")
message(STATUS "Swig: Python (${SWIG_PYTHON}) Matlab (${SWIG_MATLAB})")
message(STATUS "System name:version ${CMAKE_SYSTEM_NAME}:${CMAKE_SYSTEM_VERSION}")
message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")
message(STATUS "Installation directory is ${CMAKE_INSTALL_PREFIX}")
