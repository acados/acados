#! /bin/bash


LICENSE_TEXT_PATH="$(pwd)/utils/license/acados_license.txt"
WRONG_FOLDER_MSG="To run this script 'cd' to acados repo main directory"

[ ! -f $LICENSE_TEXT_PATH ] && echo $WRONG_FOLDER_MSG && exit

SCRIPT_REPO_URL=https://gitlab.syscop.de/tmmsartor/licenseheaders.git
SCRIPT_CMD="python licenseheaders/licenseheaders.py -D --dry"

[ ! -d licenseheaders ] && git clone $SCRIPT_REPO_URL

$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d acados/
$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d cmake/
$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d ci/
$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d interfaces/
$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d examples/
$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d test/
$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d docs/

mkdir mk

mv CMakeLists.txt Makefile Makefile.rule Makefile.osqp mk

$SCRIPT_CMD -t $LICENSE_TEXT_PATH -d mk/

pushd mk
mv CMakeLists.txt Makefile Makefile.rule Makefile.osqp ..
popd

rm -r mk
