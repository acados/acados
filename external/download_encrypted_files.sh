#!/bin/sh

INSTALL_MATLAB=0
INSTALL_OOQP_AND_DEPS=1

if [ "$INSTALL_MATLAB" = "1" ]; then
	echo "Downloading MATLAB"
	source download_matlab
else
	echo "Skipping MATLAB installation"
fi


if [ "$INSTALL_OOQP_AND_DEPS" = "1" ]; then
	echo "Downloading OOQP and dependencies"
	source encrypted_script
else
	echo "Skipping OOQP installation"
fi
