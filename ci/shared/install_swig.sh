#!/bin/bash

pushd "${TRAVIS_BUILD_DIR}/external/swig";
    ./autogen.sh;
    ./configure --prefix=$(pwd)/swig_install --enable-silent-rules;
    make;
    make install > /dev/null;    # quiet installation
    export PATH=$(pwd):$PATH;    # add swig to PATH
popd;