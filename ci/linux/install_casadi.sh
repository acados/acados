#!/bin/bash
CASADI_VERSION='3.4.0';
_CASADI_GITHUB_RELEASES="https://github.com/casadi/casadi/releases/download/${CASADI_VERSION}";
CASADI_PYTHON_URL="${_CASADI_GITHUB_RELEASES}/casadi-linux-py35-v${CASADI_VERSION}-64bit.tar.gz";
CASADI_MATLAB_URL="${_CASADI_GITHUB_RELEASES}/casadi-linux-matlabR2014b-v${CASADI_VERSION}.tar.gz";

pushd external;
	wget -O casadi-linux-py35.tar.gz "${CASADI_PYTHON_URL}";
	mkdir -p casadi-linux-py35;
	tar -xf casadi-linux-py35.tar.gz -C casadi-linux-py35;
	export PYTHONPATH=$(pwd)/casadi-linux-py35:$PYTHONPATH;
	export CASADIPATH=$(pwd)/casadi-linux-py35;

	wget -O casadi-linux-matlabR2014b.tar.gz "${CASADI_MATLAB_URL}";
	mkdir -p casadi-linux-matlabR2014b;
	tar -xf casadi-linux-matlabR2014b.tar.gz -C casadi-linux-matlabR2014b;
	export MATLABPATH=$(pwd)/casadi-linux-matlabR2014b:$MATLABPATH;
popd;