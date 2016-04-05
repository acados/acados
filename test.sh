#!/usr/bin/env bash

# fail on any error
set -e

# test hpmpc
pushd hpmpc
# TODO: hpmpc needs a test target
# make test
popd

# test acados
pushd test_problems
./test.out
popd
