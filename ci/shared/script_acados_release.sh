#!/bin/bash -e
ACADOS_INSTALL_DIR="${ACADOS_INSTALL_DIR:-${HOME}/acados}";

function build_acados {
	BUILD_TYPE="Debug";  # Release or Debug
	UNIT_TESTS='ON';
	if [ "${1}" = 'Release' ]; then
		BUILD_TYPE='Release';
		UNIT_TESTS='OFF';
	fi

	cmake -E make_directory build;
	cmake -E chdir build cmake \
			-D CMAKE_BUILD_TYPE="${BUILD_TYPE}" \
			-D UNIT_TESTS="${UNIT_TESTS}" \
			-D CMAKE_INSTALL_PREFIX="${ACADOS_INSTALL_DIR}" \
			-D SWIG_PYTHON=ON \
			-D BUILD_SHARED_LIBS=ON \
			..;
	cmake --build build;
	cmake -E chdir build ctest --output-on-failure;
	# only install release version
	[ "${BUILD_TYPE}" = 'Release' ] && cmake --build build --target install;
}

build_acados Debug;
rm -r build;
build_acados Release;