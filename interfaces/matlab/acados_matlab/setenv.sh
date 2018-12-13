# if model folder not specified assume parent of this folder
ACADOS_INSTALL_DIR=${ACADOS_INSTALL_DIR:-"$(pwd)/../../.."}

# export casadi folder
export MATLABPATH=$MATLABPATH:$ACADOS_INSTALL_DIR/external/casadi-matlabR2014b-v3.4.0/

# if model folder not specified assume this folder
MODEL_FOLDER=${MODEL_FOLDER:-"."}

# export acados and model folder
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ACADOS_INSTALL_DIR/lib:$MODEL_FOLDER
