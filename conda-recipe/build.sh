#!/bin/bash
export WITH_OPENMP=0
export BLOCK_DIAG_ILU_WITH_OPENMP=0
python -m pip install --no-deps --ignore-installed .
