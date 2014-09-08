#!/bin/bash -x
#
# Test script for the Continuous Integration service at drone.io

# Define the version of Python that should be tested
PY_VERSION="2.7"
MINICONDA_PATH=$HOME/miniconda
# Install prerequisities

# miniconda
bash scripts/install_miniconda.sh $PY_VERSION $MINICONDA_PATH
if [ $? -ne 0 ]; then
    echo "install_miniconda.sh failed"
    exit 1
fi
export PATH="$MINICONDA_PATH/bin:$PATH"

# test-env
bash scripts/conda_create_test-env.sh $PY_VERSION
if [ $? -ne 0 ]; then
    echo "conda_create_test-env failed"
    exit 1
fi

# apt-get
scripts/aptget.ubuntu.12.04LTS.sh
if [[ $? != 0 ]]; then
    echo "aptget.ubuntu.12.04LTS.sh failed."
    exit 1
fi

# sundials
bash scripts/install_sundials_w_lapack.sh
if [ $? -ne 0 ]; then
    echo "install_sundials_w_lapack.sh failed"
    exit 1
fi


# Build extension module and run test suite
source activate test-env

python -c "import sympy; print(sympy.__version__)"
python -c "import pycompilation; print(pycompilation.__version__)"
python -c "import pycodeexport; print(pycodeexport.__version__)"

export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
DISTUTILS_DEBUG=1 USE_OPENMP=1 python setup.py build_ext -i
PYTHONPATH=.:$PYTHONPATH py.test --slow --pep8 --doctest-modules --ignore build/ --ignore setup.py
if [[ $? != 0 ]]; then
    echo "py.test failed."
    exit 1
fi

CHEMREAC_SOLVER=sundials PYTHONPATH=.:$PYTHONPATH py.test --slow --ignore build/
if [[ $? != 0 ]]; then
    echo "(sundials) py.test failed."
    exit 1
fi

# Build docs
cd docs/
PYTHONPATH=`pwd`/.. make html >_build.log
cp _build.log _build/html/
cd _build/
tar -C html/ -jcf html.tar.bz2 .
