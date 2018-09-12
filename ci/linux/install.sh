#!/bin/bash -e

source "${SCRIPT_DIR}/install_apt_dependencies.sh";
source "${SCRIPT_DIR}/install_python_dependencies.sh";
source "${SHARED_SCRIPT_DIR}/install_eigen.sh";
source "${SCRIPT_DIR}/install_casadi.sh";
source "${SHARED_SCRIPT_DIR}/install_swig.sh";