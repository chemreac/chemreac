#!/bin/bash
export LIBRARY_PATH="$PREFIX/lib:$LIBRARY_PATH"
export CPLUS_INCLUDE_PATH="$PREFIX/include:$CPLUS_INCLUDE_PATH"

export WITH_OPENMP=1
export BLOCK_DIAG_ILU_WITH_DGETRF=1
export BLOCK_DIAG_ILU_WITH_OPENMP=1
python -m pip install --no-deps --ignore-installed .
