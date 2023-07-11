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

EXTERNAL_SOFTWARE_DOWNLOAD_DIR="${TRAVIS_BUILD_DIR}/external";
MATLAB_DOWNLOAD_LINK="${MATLAB_DOWNLOAD_LINK_LINUX}";
[ "${TRAVIS_OS_NAME}" = 'osx' ] && MATLAB_DOWNLOAD_LINK="${MATLAB_DOWNLOAD_LINK_OSX}";

# initialize
export MATLAB_EXECUTABLE='testing'

if [ -n "${MATLAB_DOWNLOAD_LINK}" -a "${ACADOS_MATLAB}" = 'ON' ]; then
    pushd "${EXTERNAL_SOFTWARE_DOWNLOAD_DIR}";
        wget -q -O matlab.tar.gz "${MATLAB_DOWNLOAD_LINK}";
        tar -xzf matlab.tar.gz > /dev/null
        export MATLAB_ROOT=$(pwd)/matlab
        export MATLAB_EXECUTABLE=$(pwd)/matlab/bin/matlab
    popd;
fi
