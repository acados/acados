#!/bin/bash
#
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
#
COVERAGE="${COVERAGE:-}";
ACADOS_ROOT_FOLDER="init";

export MATLABPATH="${ACADOS_INSTALL_DIR}/lib:${MATLABPATH}";

function build_acados {
	BUILD_TYPE="Debug";  # Release or Debug
	[ -z "$_JAIL" ] && echo "Empty: Yes" || echo "Empty: No"
	if [ "${1}" = 'Release' ]; then
		BUILD_TYPE='Release';
		ACADOS_LINT='OFF';
	fi

	if [ "${ACADOS_UNIT_TESTS}" = 'ON' ]; then
		ACADOS_WITH_QPOASES='ON';
	fi

	ACADOS_ROOT_FOLDER="$(pwd)";
	echo
	echo "ACADOS_ROOT_FOLDER=$ACADOS_ROOT_FOLDER" #/home/travis/build/acados/acados

	[ -d ./build ] && rm -r build;
	cmake -E make_directory build;
	cmake -E chdir build cmake \
		-D BLASFEO_TARGET="${BLASFEO_TARGET}" \
		-D HPIPM_TARGET="${HPIPM_TARGET}" \
		-D CMAKE_BUILD_TYPE="${BUILD_TYPE}" \
		-D ACADOS_UNIT_TESTS="${ACADOS_UNIT_TESTS}" \
		-D ACADOS_WITH_QPOASES="${ACADOS_WITH_QPOASES}" \
		-D ACADOS_LINT="${ACADOS_LINT}" \
		-D ACADOS_INSTALL_DIR="${ACADOS_INSTALL_DIR}" \
		-D Matlab_ROOT_DIR="${MATLAB_ROOT}" \
		-D SWIG_MATLAB="${SWIG_MATLAB}" \
		-D COVERAGE="${COVERAGE}" \
		-D SWIG_PYTHON="${SWIG_PYTHON}" \
		-D BUILD_SHARED_LIBS=ON \
		-D ACADOS_EXAMPLES="${ACADOS_EXAMPLES}" \
		-D MATLAB_EXECUTABLE="${MATLAB_EXECUTABLE}" \
		-D ACADOS_MATLAB="${ACADOS_MATLAB}" \
		-D ACADOS_OCTAVE="${ACADOS_OCTAVE}" \
		..;
	if [ "${ACADOS_LINT}" = 'ON' ]; then
		cmake --build build --target lint;
		[ $? -ne 0 ] && exit 110;
	fi


	cmake --build build;
	cmake --build build --target install;

	# echo "searching the libs"
	# find $(pwd) -name 'libhpipm.*';

	# Prepare ctest with Matlab/Octave interface
	if [[ "${ACADOS_OCTAVE}" = 'ON' || "${ACADOS_MATLAB}" = 'ON' ]]; then

		# mkdir -p /home/travis/build/acados/acados/lib;
		# echo "creating directory ${ACADOS_ROOT_FOLDER}/lib";
		# mkdir -p "${ACADOS_ROOT_FOLDER}/lib";

		# echo "creating symboic links to libaries"
		# echo
		# ln -s "${ACADOS_ROOT_FOLDER}/build/external/hpipm/libhpipm.so" "${ACADOS_ROOT_FOLDER}/lib";
		# ln -s "${ACADOS_ROOT_FOLDER}/build/external/blasfeo/libblasfeo.so" "${ACADOS_ROOT_FOLDER}/lib";
		# ln -s "${ACADOS_ROOT_FOLDER}/build/acados/libacados.so" "${ACADOS_ROOT_FOLDER}/lib";

		# Export paths
		export OCTAVE_PATH="${ACADOS_ROOT_FOLDER}/interfaces/acados_matlab":$OCTAVE_PATH;
		export ACADOS_INSTALL_DIR="${ACADOS_ROOT_FOLDER}";

		pushd examples/matlab_mex/pendulum_on_cart_model;
			# source env_ci.sh;
			MODEL_FOLDER=${MODEL_FOLDER:-"./build"}
			export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ACADOS_ROOT_FOLDER/lib:$MODEL_FOLDER
			echo
			echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
		popd;

		echo "OCTAVE_PATH=$OCTAVE_PATH";
	fi

	# Run ctest
	# TODO: test matlab/python
	cmake -E chdir build ctest --output-on-failure; # use -V for full output

	[ $? -ne 0 ] && exit 100;
	if [ -n "${COVERAGE}" ]; then
		cmake --build build --target acados_coverage || \
		  echo "Coverage report not generated";
	fi
}

# build_acados Debug;
build_acados Release;
