#! /bin/bash

SCRIPT_REPO_URL=https://gitlab.syscop.de/tmmsartor/licenseheaders.git
LICENSE_TEXT_PATH="./acados_license.txt"
SCRIPT_CMD="python licenseheaders/licenseheaders.py"

[ ! -d licenseheaders ] && git clone $SCRIPT_REPO_URL


$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d acados/
$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d cmake/
$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d ci/
$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d interfaces/
$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d examples/
$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d test/
$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d docs/

mv CMakeLists.txt Makefile Makefile.rule Makefile.osqp mk

$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d mk/

pushd mk
mv CMakeLists.txt Makefile Makefile.rule Makefile.osqp ..
popd
