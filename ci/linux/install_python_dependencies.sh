#!/bin/bash -e
virtualenv --python=python3.5 acadosenv;
source acadosenv/bin/activate;
pip install numpy scipy matplotlib;