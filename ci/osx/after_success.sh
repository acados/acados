#!/bin/bash -e

ACADOS_INSTALL_DIR="${ACADOS_INSTALL_DIR:-${HOME}/acados}";
DEPLOY_FOLDER="${DEPLOY_FOLDER:-${HOME}/deploy}";

pushd "${ACADOS_INSTALL_DIR}";
tar -zcf ${DEPLOY_FOLDER}/osx.tar.gz ./lib ./include;
popd;