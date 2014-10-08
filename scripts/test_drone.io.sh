#!/bin/bash -x -e
#
# Test script for the Continuous Integration service at drone.io

# Define the version of Python that should be tested
PY_VERSION="2.7"
MINICONDA_PATH=$HOME/miniconda
# Install prerequisities
function finish {
    # Tell github that the shield cache has expired
    curl -X PURGE https://camo.githubusercontent.com/7e92a4f31426afab773c4cc869f88e5c6a38e946/68747470733a2f2f64726f6e652e696f2f6769746875622e636f6d2f626a6f6461682f6368656d726561632f7374617475732e706e67
}
trap finish EXIT

# miniconda
bash scripts/install_miniconda.sh $PY_VERSION $MINICONDA_PATH
export PATH="$MINICONDA_PATH/bin:$PATH"

# test-env
bash scripts/setup_conda_testenv.sh $PY_VERSION $MINICONDA_PATH $ENV_NAME

# apt-get
scripts/ubuntu12.04/apt_get_gcc_48.sh
scripts/ubuntu12.04/apt_get_g++_48.sh
scripts/ubuntu12.04/apt_get_gfortran_48.sh
scripts/aptget_debian.sh

# sundials
bash scripts/install_sundials_w_lapack.sh

mkdir -p $HOME/.config/matplotlib/
cp ./scripts/matplotlibrc $HOME/.config/matplotlib/

# Build extension module and run test suite
source activate test-env

python -c "import sympy; print(sympy.__version__)"
python -c "import pycompilation; print(pycompilation.__version__)"
python -c "import pycodeexport; print(pycodeexport.__version__)"

export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
DISTUTILS_DEBUG=1 USE_OPENMP=1 python setup.py build_ext -i
PYTHONPATH=`pwd`:$PYTHONPATH ./scripts/run_tests.sh
if [[ $? != 0 ]]; then
    echo "run_tests.sh failed."
    exit 1
fi
tar -jcf htmlcov.tar.bz2 htmlcov/

CHEMREAC_SOLVER=sundials PYTHONPATH=.:$PYTHONPATH py.test --slow --veryslow --ignore build/

# Build docs
bash scripts/build_docs.sh
tar -jcf html.tar.bz2 docs/_build/html/

