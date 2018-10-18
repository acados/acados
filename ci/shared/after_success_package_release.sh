#!/bin/bash

ACADOS_INSTALL_DIR="${ACADOS_INSTALL_DIR:-${HOME}/acados}";
DEPLOY_FOLDER="${DEPLOY_FOLDER:-${HOME}/deploy}";
DEPLOY_NAME="${DEPLOY_NAME:-default}";

pushd "${ACADOS_INSTALL_DIR}";
tar -zcf ${DEPLOY_FOLDER}/${DEPLOY_NAME}.tar.gz ./lib ./include;
popd;