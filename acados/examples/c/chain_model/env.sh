#! /usr/bin/bash
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


#if [[ "${BASH_SOURCE[0]}" != "${0}" ]]
#then
#	echo "Script is being sourced"
#else
#	echo "ERROR: Script is a subshell"
#	echo "To affect your current shell enviroment source this script with:"
#	echo "source env.sh"
#	exit
#fi

# check that this file is run
#export ENV_RUN=true

# if acados folder not specified assume parent of the folder of the single examples
ACADOS_INSTALL_DIR=${ACADOS_INSTALL_DIR:-"$(pwd)/../../.."}
export ACADOS_INSTALL_DIR
echo
echo "ACADOS_INSTALL_DIR=$ACADOS_INSTALL_DIR"

# export casadi folder and matlab/octave mex folder
# MATLAB case
export MATLABPATH=$MATLABPATH:$ACADOS_INSTALL_DIR/external/casadi-matlab/:$ACADOS_INSTALL_DIR/interfaces/acados_matlab_octave/
echo
echo "MATLABPATH=$MATLABPATH"
# Octave case
export OCTAVE_PATH=$OCTAVE_PATH:$ACADOS_INSTALL_DIR/external/casadi-octave/:$ACADOS_INSTALL_DIR/interfaces/acados_matlab_octave/
echo
echo "OCTAVE_PATH=$OCTAVE_PATH"

# export acados mex flags
#export ACADOS_MEX_FLAGS="GCC=/usr/bin/gcc-4.9"

# if model folder not specified assume this folder
#MODEL_FOLDER=${MODEL_FOLDER:-"./build"}
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ACADOS_INSTALL_DIR/lib:$MODEL_FOLDER
#echo
#echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
