#!/bin/bash
#
# Test script for the Continuous Integration service at drone.io

# Define the version of Python that should be tested
PYTHON_VERSION="2.7"
CONDA_PY="27"
ENV_NAME="testenv"
MINICONDA_PATH=$HOME/miniconda
export DEBIAN_FRONTEND=noninteractive
# Install prerequisities
function finish {
    # Tell github that the shield cache has expired
    curl -X PURGE https://camo.githubusercontent.com/7e92a4f31426afab773c4cc869f88e5c6a38e946/68747470733a2f2f64726f6e652e696f2f6769746875622e636f6d2f626a6f6461682f6368656d726561632f7374617475732e706e67
}
trap finish EXIT

# install miniconda in the background in a subshell
(
source scripts/install_miniconda.sh $PYTHON_VERSION $MINICONDA_PATH "3.10.1"
conda config --add channels http://conda.anaconda.org/bjodah
) &
miniconda_install_pid=$!

# Install compilers and dependencies
scripts/ubuntu12.04/apt_get_gcc_48.sh
scripts/ubuntu12.04/apt_get_g++_48.sh
scripts/ubuntu12.04/apt_get_gfortran_48.sh
scripts/aptget_debian.sh

# Download and install sundials
bash scripts/install_sundials_w_lapack.sh
export LIBRARY_PATH=/usr/local/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

# Set non-interactive matplotlib backend
mkdir -p $HOME/.config/matplotlib/
cp ./scripts/matplotlibrc $HOME/.config/matplotlib/

# Wait for miniconda installation
wait $miniconda_install_pid
export PATH="$MINICONDA_PATH/bin:$PATH"
hash -r

# Build extension module and run test suite
export WITH_OPENMP=1
export CHEMREAC_INTEGRATION_KWARGS="{'integrator': cvode, 'linear_solver': 'direct'}"
source ./scripts/ci_conda.sh $PYTHON_VERSION $CONDA_PY $ENV_NAME 1
if [[ $? != 0 ]]; then
    >&2 echo "./scripts/ci_conda.sh failed."
    exit 1
fi
# Archive the coverage report (for artifact)
tar -jcf htmlcov.tar.bz2 htmlcov/

# Build docs and archive (for artifact)
bash scripts/build_docs.sh
tar -jcf html.tar.bz2 docs/_build/html/
