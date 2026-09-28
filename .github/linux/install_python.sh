#!/bin/bash
#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

python3 --version;

# install virtualenv
sudo pip3 install virtualenv;
# create virtualenv
virtualenv --python=python3 acadosenv;
# source virtualenv
source acadosenv/bin/activate;
which python;

pip install interfaces/acados_template
