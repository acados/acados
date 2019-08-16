#!/bin/bash
#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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