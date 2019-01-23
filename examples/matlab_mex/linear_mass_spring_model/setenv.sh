# if acados folder not specified assume parent of this folder
ACADOS_FOLDER=${ACADOS_FOLDER:-"$(pwd)/../../.."}
export ACADOS_FOLDER

# export casadi folder and matlab mex folder
export MATLABPATH=$MATLABPATH:$ACADOS_FOLDER/external/casadi-matlabR2014b-v3.4.0/:$ACADOS_FOLDER/interfaces/acados_matlab/

# export acados mex flags
#export ACADOS_MEX_FLAGS="GCC=/usr/bin/gcc-4.9"

# if model folder not specified assume this folder
MODEL_FOLDER=${MODEL_FOLDER:-"."}
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ACADOS_FOLDER/lib:$MODEL_FOLDER
