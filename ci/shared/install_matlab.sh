#/bin/bash
EXTERNAL_SOFTWARE_DOWNLOAD_DIR="${TRAVIS_BUILD_DIR}/external";
MATLAB_DOWNLOAD_LINK="${MATLAB_DOWNLOAD_LINK_LINUX}";
[ "${TRAVIS_OS_NAME}" = 'osx' ] && MATLAB_DOWNLOAD_LINK="${MATLAB_DOWNLOAD_LINK_OSX}";

if [ -n "${MATLAB_DOWNLOAD_LINK}" -a "${SWIG_MATLAB}" = 'ON' ]; then
    pushd "${EXTERNAL_SOFTWARE_DOWNLOAD_DIR}";
        wget -q -O matlab.tar.gz "${MATLAB_DOWNLOAD_LINK}";
        tar -xzf matlab.tar.gz > /dev/null
        export MATLAB_ROOT=$(pwd)/matlab
    popd;
fi
