#!/bin/bash

wget "http://hera.physchem.kth.se/~chemreac/$1.tar.gz"
tar xvzf "$1.tar.gz"
pushd "$1"
CPLUS_INCLUDE_PATH=/usr/include/python2.7 ./scripts/run_tests.sh
popd
