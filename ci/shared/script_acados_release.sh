#!/bin/bash
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

# COVERAGE="${COVERAGE:-}";


export MATLABPATH="${ACADOS_INSTALL_DIR}/lib:${MATLABPATH}";

function build_acados {
	BUILD_TYPE="Debug";  # Release or Debug
	[ -z "$_JAIL" ] && echo "Empty: Yes" || echo "Empty: No"
	if [ "${1}" = 'Release' ]; then
		BUILD_TYPE='Release';
		ACADOS_LINT='OFF';
	fi

	if [[ "${ACADOS_UNIT_TESTS}" = 'ON' || "${ACADOS_PYTHON}" = 'ON' ]]; then
		ACADOS_WITH_QPOASES='ON';
	fi

	[ -d ./build ] && rm -r build;
	cmake -E make_directory build;
	cmake -E chdir build cmake \
		-D BLASFEO_TARGET="${BLASFEO_TARGET}" \
		-D HPIPM_TARGET="${HPIPM_TARGET}" \
		-D CMAKE_BUILD_TYPE="${BUILD_TYPE}" \
		-D ACADOS_UNIT_TESTS="${ACADOS_UNIT_TESTS}" \
		-D ACADOS_WITH_QPOASES="${ACADOS_WITH_QPOASES}" \
		-D ACADOS_WITH_QPDUNES="${ACADOS_WITH_QPDUNES}" \
		-D ACADOS_WITH_OSQP="${ACADOS_WITH_OSQP}" \
		-D ACADOS_LINT="${ACADOS_LINT}" \
		-D ACADOS_INSTALL_DIR="${ACADOS_INSTALL_DIR}" \
		-D Matlab_ROOT_DIR="${MATLAB_ROOT}" \
		-D COVERAGE="${COVERAGE}" \
		-D BUILD_SHARED_LIBS=ON \
		-D ACADOS_EXAMPLES="${ACADOS_EXAMPLES}" \
		-D MATLAB_EXECUTABLE="${MATLAB_EXECUTABLE}" \
		-D ACADOS_MATLAB="${ACADOS_MATLAB}" \
		-D ACADOS_OCTAVE="${ACADOS_OCTAVE}" \
		-D ACADOS_PYTHON="${ACADOS_PYTHON}" \
		..;
	[ $? -ne 0 ] && exit 110;

	if [ "${ACADOS_LINT}" = 'ON' ]; then
		cmake --build build --target lint;
		[ $? -ne 0 ] && exit 110;
	fi

	cmake --build build;
	[ $? -ne 0 ] && exit 110;

	cmake --build build --target install;
	[ $? -ne 0 ] && exit 110;

    if [[ "${ACADOS_PYTHON}" = 'ON' ]] ;
    then
        source "${SCRIPT_DIR}/install_python_dependencies.sh";
        pushd interfaces/acados_template;
            pip install .
        popd;
        source "${SCRIPT_DIR}/install_t_renderer.sh";
    fi

	# Run ctest
	cmake -E chdir build ctest -V; # use -V for full output # --output-on-failure for less

	[ $? -ne 0 ] && exit 100;
	if [ -n "${COVERAGE}" ]; then
		echo "analyzing test coverage";
		cmake --build build --target acados_coverage || \
		  echo "Coverage report not generated";
	fi
}

# build_acados Debug;
build_acados Release;
