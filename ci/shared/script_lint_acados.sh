#!/bin/bash -e

cmake -E make_directory build;
cmake --build build --target lint;