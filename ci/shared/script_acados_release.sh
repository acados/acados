#!/bin/bash
ACADOS_INSTALL_DIR="${ACADOS_INSTALL_DIR:-${HOME}/acados}";
COVERAGE="${COVERAGE:-}";

export MATLABPATH="${ACADOS_INSTALL_DIR}/lib:${MATLABPATH}";

function build_acados {
	BUILD_TYPE="Debug";  # Release or Debug
	[ -z "$_JAIL" ] && echo "Empty: Yes" || echo "Empty: No"
	if [ "${1}" = 'Release' ]; then
		BUILD_TYPE='Release';
		ACADOS_LINT='OFF';
	fi

	[ -d ./build ] && rm -r build;
	cmake -E make_directory build;
	cmake -E chdir build cmake \
			-D CMAKE_BUILD_TYPE="${BUILD_TYPE}" \
			-D ACADOS_UNIT_TESTS="${ACADOS_UNIT_TESTS}" \
			-D ACADOS_LINT="${ACADOS_LINT}" \
			-D ACADOS_INSTALL_DIR="${ACADOS_INSTALL_DIR}" \
			-D Matlab_ROOT_DIR="${MATLAB_ROOT}" \
			-D SWIG_MATLAB="${SWIG_MATLAB}" \
			-D COVERAGE="${COVERAGE}" \
			-D SWIG_PYTHON="${SWIG_PYTHON}" \
			-D BUILD_SHARED_LIBS=ON \
			..;
	if [ "${ACADOS_LINT}" = 'ON' ]; then
		cmake --build build --target lint;
		[ $? -ne 0 ] && exit 110;
	fi

	cmake --build build;
	cmake -E chdir build ctest --output-on-failure;
	[ $? -ne 0 ] && exit 100;
	if [ -n "${COVERAGE}" ]; then
		cmake --build build --target acados_coverage || \
		  echo "Coverage report not generated";
	fi
	# only install release version
	[ "${BUILD_TYPE}" = 'Release' ] && cmake --build build --target install;
}

build_acados Debug;
build_acados Release;
