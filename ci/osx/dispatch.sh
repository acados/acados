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


if [ "${SECTION}" = 'before_install' ]; then
    # export GENERIC target for osx
    export BLASFEO_TARGET=GENERIC;
    export HPIPM_TARGET=GENERIC;

    export ACADOS_INSTALL_DIR="$(pwd)";
    export ACADOS_SOURCE_DIR="$(pwd)";

elif [ "${SECTION}" = 'install' ]; then
    source "${SCRIPT_DIR}/install_ccache.sh";
    source "${SHARED_SCRIPT_DIR}/install_eigen.sh";
    source "${SCRIPT_DIR}/install_python.sh";

    if  [[ "${ACADOS_MATLAB}" = 'ON' || "${ACADOS_OCTAVE}" = 'ON' ]] ||
        [[ "${ACADOS_PYTHON}" = 'ON' ]];
        then
        source "${SCRIPT_DIR}/install_casadi.sh";
    fi

    if [[ "${ACADOS_PYTHON}" = 'ON' ]] ;
    then
        source "${SCRIPT_DIR}/install_python_dependencies.sh";
        pushd examples/acados_template/python/test;
            export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ACADOS_INSTALL_DIR/lib:$MODEL_FOLDER
        popd;
    fi

    if [[ "${ACADOS_MATLAB}" = 'ON' ]];
    then
        source "${SHARED_SCRIPT_DIR}/install_matlab.sh";
    fi

elif [ "${SECTION}" = 'script' ]; then
    source "${SHARED_SCRIPT_DIR}/script_acados_release.sh";

elif [ "${SECTION}" = 'after_success' ]; then
    source "${SHARED_SCRIPT_DIR}/after_success_package_release.sh";

fi
