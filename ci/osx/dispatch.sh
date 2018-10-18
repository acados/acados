#!/bin/bash

if [ "${SECTION}" = 'install' ]; then
	source "${SCRIPT_DIR}/install_ccache.sh";
	source "${SCRIPT_DIR}/install_python_dependencies.sh";
	source "${SHARED_SCRIPT_DIR}/install_eigen.sh";
	source "${SHARED_SCRIPT_DIR}/install_matlab.sh"
	source "${SHARED_SCRIPT_DIR}/install_swig.sh";
	source "${SCRIPT_DIR}/install_casadi.sh";

elif [ "${SECTION}" = 'script' ]; then
	source "${SHARED_SCRIPT_DIR}/script_acados_release.sh";

elif [ "${SECTION}" = 'after_success' ]; then
	source "${SHARED_SCRIPT_DIR}/after_success_package_release.sh";

fi
