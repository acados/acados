#!/bin/bash -e
ACADOS_INSTALL_DIR="${ACADOS_INSTALL_DIR:-${HOME}/acados}";

cmake -E make_directory build;
cmake -E chdir build cmake \
		-D CMAKE_BUILD_TYPE=Release \
		-D CMAKE_INSTALL_PREFIX="${ACADOS_INSTALL_DIR}" \
		-D SWIG_PYTHON=ON \
		..;
cmake --build build --target install;