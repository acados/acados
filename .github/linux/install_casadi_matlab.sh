#!/bin/bash
#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

# CASADI_VERSION='3.5.3'; # Latest version with Octave 4.2.2 binaries
CASADI_VERSION='3.8.1';
MATLAB_VERSION='2018b';

_CASADI_GITHUB_RELEASES="https://github.com/casadi/casadi/releases/download/${CASADI_VERSION}";

CASADI_MATLAB_URL="${_CASADI_GITHUB_RELEASES}/casadi-${CASADI_VERSION}-linux64-matlab${MATLAB_VERSION}.zip";

wget -O casadi-linux-matlab.zip "${CASADI_MATLAB_URL}";
mkdir -p casadi-matlab;
unzip casadi-linux-matlab.zip -d ./casadi-matlab;
