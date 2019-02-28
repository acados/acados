#! /usr/bin/bash

if [[ "${BASH_SOURCE[0]}" != "${0}" ]]
then
	echo "Script is being sourced"
else
	echo "ERROR: Script is a subshell"
	echo "To affect your current shell enviroment source this script with:"
	echo "source env.sh"
	exit
fi

# if acados folder not specified assume parent of this folder
ACADOS_INSTALL_DIR=${ACADOS_INSTALL_DIR:-"$(pwd)/../../.."}
export ACADOS_INSTALL_DIR
echo
echo "ACADOS_INSTALL_DIR=$ACADOS_INSTALL_DIR"

# export casadi folder and matlab/octave mex folder
# MATLAB case
export MATLABPATH=$MATLABPATH:$ACADOS_INSTALL_DIR/external/casadi-matlabR2014b-v3.4.0/:$ACADOS_INSTALL_DIR/interfaces/acados_matlab/
echo
echo "MATLABPATH=$MATLABPATH"
# Octave case
export OCTAVE_PATH=$OCTAVE_PATH:$ACADOS_INSTALL_DIR/external/casadi-octave-v3.4.0/:$ACADOS_INSTALL_DIR/interfaces/acados_matlab/
echo
echo "OCTAVE_PATH=$OCTAVE_PATH"

# export acados mex flags
#export ACADOS_MEX_FLAGS="GCC=/usr/bin/gcc-4.9"

# if model folder not specified assume this folder
MODEL_FOLDER=${MODEL_FOLDER:-"."}
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ACADOS_INSTALL_DIR/lib:$MODEL_FOLDER
echo
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
