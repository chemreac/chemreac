#!/bin/bash
#
# Test script for the Continuous Integration service at drone.io

# Define the version of Python that should be tested
PYTHON_VERSION="2.7"
MINICONDA_PATH=$HOME/miniconda
# Install prerequisities
function finish {
    # Tell github that the shield cache has expired
    curl -X PURGE https://camo.githubusercontent.com/7e92a4f31426afab773c4cc869f88e5c6a38e946/68747470733a2f2f64726f6e652e696f2f6769746875622e636f6d2f626a6f6461682f6368656d726561632f7374617475732e706e67
}
trap finish EXIT

# miniconda
source scripts/install_miniconda.sh $PYTHON_VERSION $MINICONDA_PATH "3.7.0"
conda config --add channels http://conda.binstar.org/$BINSTAR_USER

# apt-get
scripts/ubuntu12.04/apt_get_gcc_48.sh
scripts/ubuntu12.04/apt_get_g++_48.sh
scripts/ubuntu12.04/apt_get_gfortran_48.sh
scripts/aptget_debian.sh

# sundials
bash scripts/install_sundials_w_lapack.sh
export LIBRARY_PATH=/usr/local/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

mkdir -p $HOME/.config/matplotlib/
cp ./scripts/matplotlibrc $HOME/.config/matplotlib/

# Build extension module and run test suite
export USE_OPENMP=1
source ./scripts/ci_conda.sh $PYTHON_VERSION $CONDA_PY $ENV_NAME 1
if [[ $? != 0 ]]; then
    >&2 echo "./scripts/ci_conda.sh failed."
    exit 1
fi
tar -jcf htmlcov.tar.bz2 htmlcov/

CHEMREAC_SOLVER=sundials py.test --slow --veryslow --ignore build/

# Build docs
bash scripts/build_docs.sh
tar -jcf html.tar.bz2 docs/_build/html/
