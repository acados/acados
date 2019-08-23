#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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


set(CMAKE_SYSTEM_NAME dSpace)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

file(TO_CMAKE_PATH "C:\\Program Files\\dSPACE RCPHIL 2016-B" DSPACE_TOOLS)
set(DSPACE_RTLIB "${DSPACE_TOOLS}/DS1401/RTLib")
set(DSPACE_PPCTOOLS ${DSPACE_TOOLS}/Compiler/PPCTools)
find_program(CMAKE_C_COMPILER NAMES ${DSPACE_PPCTOOLS}/bin/mccppc.exe)
find_program(CMAKE_CXX_COMPILER NAMES ${DSPACE_PPCTOOLS}/bin/cccppc.exe)
find_program(CMAKE_ASM_COMPILER NAMES ${DSPACE_PPCTOOLS}/bin/asmppc.exe)
find_program(CMAKE_AR NAMES ${DSPACE_PPCTOOLS}/bin/cccppc.exe)
set(CMAKE_RANLIB ":")

set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> <OBJECTS> -o <TARGET>")
set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> <OBJECTS> -o <TARGET>")

find_program(CMAKE_MAKE_PROGRAM NAMES ${DSPACE_TOOLS}/Exe/DSMAKE.EXE)

set(CMAKE_FIND_ROOT_PATH ${DSPACE_PPCTOOLS})

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
