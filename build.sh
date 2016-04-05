#!/usr/bin/env bash

# fail on any error
set -e

# build hpmpc
pushd hpmpc
make
popd

# build test problems
pushd test_problems
make
popd
