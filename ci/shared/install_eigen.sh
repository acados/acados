#!/bin/bash
EIGEN_URL="http://bitbucket.org/eigen/eigen/get/3.2.10.tar.gz";

pushd external;
	wget -O eigen.tar.gz "${EIGEN_URL}";
	mkdir -p eigen;
	tar -xf eigen.tar.gz --strip-components=1 -C eigen;
popd;