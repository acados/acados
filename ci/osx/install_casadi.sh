#!/bin/bash
CASADIPATH="${CASADIPATH:-${HOME}/casadi}";
CASADI_VERSION='3.4.0';
CASADI_DOWNLOAD_URL="https://github.com/casadi/casadi/archive/${CASADI_VERSION}.zip";

# run only if casadi build was not cached
if [ ! -d "${CASADIPATH}" -o -z "$(ls -A "${CASADIPATH}")" ]; then
	pushd "${TRAVIS_BUILD_DIR}/external";
		curl -o casadi.zip -Ls "${CASADI_DOWNLOAD_URL}";
		unzip -qq casadi.zip;
		rm -f casadi.zip;
		pushd "./casadi-${CASADI_VERSION}";
			mkdir build;
			pushd ./build;
				cmake -DWITH_SELFCONTAINED=ON -DWITH_PYTHON=ON -DWITH_PYTHON3=ON -DCMAKE_INSTALL_PREFIX="${CASADIPATH}" ..;
				make -j 4 && make install;
			popd;
		popd;
	popd;
fi
export PYTHONPATH="${CASADIPATH}:$PYTHONPATH";