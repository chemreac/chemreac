#!/bin/bash

wget "http://hera.physchem.kth.se/~bjorn/chemreac/$1.tar.gz"
tar xvzf "$1.tar.gz"
pushd "$1"
#pip install -r requirements.txt
./scripts/run_tests.sh
popd
