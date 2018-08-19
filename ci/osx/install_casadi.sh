#!/bin/bash -e
CASADI_INSTALL_DIR="${CASADI_INSTALL_DIR}";
CASADI_VERSION='3.4.0';
CASADI_DOWNLOAD_URL="https://github.com/casadi/casadi/archive/${CASADI_VERSION}.zip";

# run only if casadi build was not cached
if [ ! -d "${CASADI_INSTALL_DIR}" ]; then
	pushd "${TRAVIS_BUILD_DIR}/external";
		curl -o casadi.zip -Ls "${CASADI_DOWNLOAD_URL}";
		unzip -qq casadi.zip;
		rm -f casadi.zip;
		pushd "./casadi-${CASADI_VERSION}";
			mkdir build;
			pushd ./build;
				cmake -DWITH_SELFCONTAINED=ON -DWITH_PYTHON=ON -DWITH_PYTHON3=ON -DCMAKE_INSTALL_PREFIX="${CASADI_INSTALL_DIR}" ..;
				make -j 4 && make install;
			popd;
		popd;
	popd;
fi
export PYTHONPATH="${CASADI_INSTALL_DIR}:$PYTHONPATH";
export CASADIPATH="${CASADI_INSTALL_DIR}";