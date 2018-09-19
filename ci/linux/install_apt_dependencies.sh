#!/bin/bash
sudo apt-get update -yqq;
sudo apt-get --allow-unauthenticated install -yqq $CXX $CC $COVERAGE python3.5 python3.5-tk python3.5-dev;