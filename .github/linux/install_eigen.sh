#!/bin/bash
#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

EIGEN_URL="https://gitlab.com/libeigen/eigen/-/archive/3.2.10/eigen-3.2.10.tar.gz";

pushd external;
	wget -O eigen.tar.gz "${EIGEN_URL}";
	mkdir -p eigen;
	tar -xf eigen.tar.gz --strip-components=1 -C eigen;
popd;