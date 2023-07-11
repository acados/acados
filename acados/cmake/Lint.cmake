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


# Lint
find_package(PythonInterp 3)

if(NOT PYTHONINTERP_FOUND)
    message(STATUS "Python Interpreter not found, disabling lint")
    return()
endif()

set(FIND_FILES_TO_LINT
	find acados examples interfaces test
	-type f -name "*.c" -o -name "*.cpp" -o -name "*.h" -o -name "*.hpp" -o -name "*.i")
set(FIND_FILES_TO_LINT ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${FIND_FILES_TO_LINT})

execute_process(COMMAND ${FIND_FILES_TO_LINT} OUTPUT_VARIABLE FILES_TO_LINT)
string(REPLACE "\n" " " FILES_TO_LINT ${FILES_TO_LINT})
separate_arguments(FILES_TO_LINT)

set(LINT_COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}
	${PYTHON_EXECUTABLE} ./utils/cpplint.py --quiet --counting=detailed
	--extensions=c,cpp,h,hpp,i ${FILES_TO_LINT})

add_custom_target(lint ${LINT_COMMAND})
