#!/bin/bash
#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

# CASADI_VERSION='3.5.3'; # Latest version with Octave 4.2.2 binaries
CASADI_VERSION='3.8.1';
OCTAVE_VERSION='7.3.0';

_CASADI_GITHUB_RELEASES="https://github.com/casadi/casadi/releases/download/${CASADI_VERSION}";

CASADI_OCTAVE_URL="${_CASADI_GITHUB_RELEASES}/casadi-${CASADI_VERSION}-linux64-octave${OCTAVE_VERSION}.zip";

# URL for Octave new CasADi
# CASADI_OCTAVE_URL="https://github.com/casadi/casadi/releases/download/nightly-se/casadi-se-linux64-octave7.3.0.zip"

wget -O casadi-linux-octave.zip "${CASADI_OCTAVE_URL}";
mkdir -p casadi-octave;
unzip casadi-linux-octave.zip -d ./casadi-octave;
