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

echo "ACADOS_SOURCE_DIR=$1/acados" >> $GITHUB_ENV
echo "ACADOS_INSTALL_DIR=$1/acados" >> $GITHUB_ENV
echo "LD_LIBRARY_PATH=$1/acados/lib" >> $GITHUB_ENV
echo "MATLABPATH=$MATLABPATH:$1/acados/interfaces/acados_matlab_octave:$1/acados/interfaces/acados_matlab_octave/acados_template_mex:${1}/acados/external/casadi-matlab" >> $GITHUB_ENV
echo "OCTAVE_PATH=$OCTAVE_PATH:${1}/acados/interfaces/acados_matlab_octave:${1}/acados/interfaces/acados_matlab_octave/acados_template_mex:${1}/acados/external/casadi-octave" >> $GITHUB_ENV
echo "LD_RUN_PATH=${1}/acados/examples/acados_matlab_octave/test/c_generated_code:${1}/acados/examples/acados_matlab_octave/pendulum_on_cart_model/c_generated_code:${1}/acados/examples/acados_matlab_octave/getting_started/c_generated_code:${1}/acados/examples/acados_matlab_octave/mocp_transition_example/c_generated_code:${1}/acados/examples/acados_matlab_octave/simple_dae_model/c_generated_code:${1}/acados/examples/acados_matlab_octave/lorentz/c_generated_code:${1}/acados/examples/acados_python/p_global_example/c_generated_code:${1}/acados/examples/acados_python/p_global_example/c_generated_code_single_phase:${1}/acados/examples/acados_python/pendulum_on_cart/sim/c_generated_code" >> $GITHUB_ENV
echo "ENV_RUN=true" >> $GITHUB_ENV
