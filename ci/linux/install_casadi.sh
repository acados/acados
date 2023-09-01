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

CASADI_VERSION='3.5.5';
OCTAVE_VERSION='5.2.0';
MATLAB_VERSION='R2014b';
_CASADI_GITHUB_RELEASES="https://github.com/casadi/casadi/releases/download/${CASADI_VERSION}";
CASADI_PYTHON_URL="${_CASADI_GITHUB_RELEASES}/casadi-linux-py35-v${CASADI_VERSION}-64bit.tar.gz";
CASADI_MATLAB_URL="${_CASADI_GITHUB_RELEASES}/casadi-linux-matlab${MATLAB_VERSION}-v${CASADI_VERSION}.tar.gz";
CASADI_OCTAVE_URL="${_CASADI_GITHUB_RELEASES}/casadi-linux-octave-${OCTAVE_VERSION}-v${CASADI_VERSION}.tar.gz";

echo "installing CasADi"

pushd external;
	if [[ "${ACADOS_PYTHON}" = 'ON' ]] ;
	then
		wget -O casadi-linux-py35.tar.gz "${CASADI_PYTHON_URL}";
		mkdir -p casadi-linux-py35;
		tar -xf casadi-linux-py35.tar.gz -C casadi-linux-py35;
		export PYTHONPATH=$(pwd)/casadi-linux-py35:$PYTHONPATH;
	fi

	if [[ "${ACADOS_MATLAB}" = 'ON' ]];
	then
		wget -O casadi-linux-matlabR2014b.tar.gz "${CASADI_MATLAB_URL}";
		mkdir -p casadi-linux-matlabR2014b;
		tar -xf casadi-linux-matlabR2014b.tar.gz -C casadi-linux-matlabR2014b;
		export MATLABPATH=$(pwd)/casadi-linux-matlabR2014b:$MATLABPATH;
	fi

	if [["${ACADOS_OCTAVE}" = 'ON' ]];
	then
		wget -O casadi-linux-octave.tar.gz "${CASADI_OCTAVE_URL}";
		mkdir -p casadi-octave;
		tar -xf casadi-linux-octave.tar.gz -C casadi-octave;
		export OCTAVE_PATH=$(pwd)/casadi-octave:$OCTAVE_PATH;
	fi
popd;
