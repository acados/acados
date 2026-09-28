#!/bin/bash
#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

echo "ACADOS_SOURCE_DIR=$1" >> $GITHUB_ENV
echo "ACADOS_INSTALL_DIR=$1" >> $GITHUB_ENV
echo "LD_LIBRARY_PATH=$1/lib:$LD_LIBRARY_PATH" >> $GITHUB_ENV
echo "MATLABPATH=$MATLABPATH:$1/interfaces/acados_matlab_octave:$1//interfaces/acados_matlab_octave/acados_template_mex:${1}/external/casadi-matlab" >> $GITHUB_ENV
echo "OCTAVE_PATH=$OCTAVE_PATH:${1}/interfaces/acados_matlab_octave:${1}/interfaces/acados_matlab_octave/acados_template_mex:${1}/external/casadi-octave" >> $GITHUB_ENV
echo "LD_RUN_PATH=${1}/examples/acados_matlab_octave/test/c_generated_code:${1}/examples/acados_matlab_octave/pendulum_on_cart_model/c_generated_code:${1}/examples/acados_matlab_octave/getting_started/c_generated_code:${1}/examples/acados_matlab_octave/mocp_transition_example/c_generated_code:${1}/examples/acados_matlab_octave/p_global_example/c_generated_code:${1}/examples/acados_matlab_octave/simple_dae_model/c_generated_code:${1}/examples/acados_matlab_octave/lorentz/c_generated_code:${1}/examples/acados_python/p_global_example/c_generated_code:${1}/examples/acados_python/p_global_example/c_generated_code_single_phase:${1}/examples/acados_python/pendulum_on_cart/sim/c_generated_code" >> $GITHUB_ENV
echo "ENV_RUN=true" >> $GITHUB_ENV
