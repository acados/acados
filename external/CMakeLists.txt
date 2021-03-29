#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren, Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor, Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan, Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-error -w")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error -w")
endif()

if(CMAKE_SYSTEM_NAME MATCHES "dSpace")
    set(USE_C99_MATH OFF CACHE BOOL "fmax is not implemented for dSpace.")
endif()

# When BLASFEO is used as part of acados
set(BLASFEO_HEADERS_INSTALLATION_DIRECTORY "include/blasfeo/include" CACHE STRING "")

set(TARGET ${BLASFEO_TARGET} CACHE STRING "Set CPU architecture target" FORCE)
set(MF "PANELMAJ" CACHE STRING "Matrix format" FORCE)
set(BLAS_API OFF CACHE BOOL "Compile BLAS API" FORCE)
add_subdirectory(blasfeo)

# When HPIPM is used as part of acados
set(HPIPM_BLASFEO_PATH ${BLASFEO_SRC_DIR} CACHE STRING "Set CPU architecture target" FORCE)

set(HPIPM_BLASFEO_LIB OFF CACHE BOOL "Turn off BLASFEO from HPIPM" FORCE)
set(HPIPM_HEADERS_INSTALLATION_DIRECTORY "include/hpipm/include" CACHE STRING "")
set(TARGET ${HPIPM_TARGET} CACHE STRING "Set CPU architecture target" FORCE)
set(HPIPM_TESTING OFF CACHE BOOL "Tests disabled" FORCE)
add_subdirectory(hpipm)

if(ACADOS_WITH_HPMPC)
	# TODO: The following definition should go ASAP!
	# When hpmpc is used as part of acados
	set(HPMPC_HEADERS_INSTALLATION_DIRECTORY "include/hpmpc/include" CACHE STRING "")
	add_subdirectory(hpmpc)
endif()

if(ACADOS_WITH_QORE)
    set(QORE_DENSE TRUE CACHE BOOL "Build dense qore library")
    add_subdirectory(qore)
endif()

if(ACADOS_WITH_QPDUNES)
    add_subdirectory(qpdunes)
endif()

if(ACADOS_WITH_QPOASES)
    add_subdirectory(qpoases)
endif()

if(ACADOS_WITH_OSQP)
    set(PROFILING OFF CACHE BOOL "Turn off profiling in OSQP")
    set(CTRLC OFF CACHE BOOL "Turn off CTRLC in OSQP")
    add_subdirectory(osqp)
endif()

if (ACADOS_WITH_OOQP)
    include(external/ooqp)
endif()
