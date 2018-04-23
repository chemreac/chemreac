#!/bin/bash
export LIBRARY_PATH="$PREFIX/lib"
export CPLUS_INCLUDE_PATH="$PREFIX/include"

export WITH_OPENMP=1
export BLOCK_DIAG_ILU_WITH_DGETRF=1
export BLOCK_DIAG_ILU_WITH_OPENMP=1
export CC=clang-6.0
export CXX=clang++
export CFLAGS=-stdlib=libc++
${PYTHON} -m pip install --no-deps --ignore-installed .
