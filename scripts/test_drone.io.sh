#!/bin/bash -x
#
# Test script for the Continuous Integration service at drone.io

# Define the version of Python that should be tested
PY_VERSION="2.7"

scripts/aptget.ubuntu.12.04LTS.sh
bash scripts/install_sundials_w_lapack.sh
if [ $? -ne 0 ]; then
    exit 1
fi

if [[ "$PY_VERSION" == "2.7" ]]; then
    # This saves us some downloading for this version
    wget --quiet http://repo.continuum.io/miniconda/Miniconda-3.6.0-Linux-x86_64.sh -O miniconda.sh;
else
    wget --quiet http://repo.continuum.io/miniconda/Miniconda3-3.6.0-Linux-x86_64.sh -O miniconda.sh;
fi
if [[ $? != 0 ]]; then
    echo "Failed to get Miniconda."
    exit 1
fi

chmod +x miniconda.sh
./miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update --quiet conda
conda info -a
conda create --quiet -n test-env python=${PY_VERSION} numpy=1.8.1 scipy=0.14 cython=0.20.1 pip
source activate test-env
pip install --quiet argh mako quantities pytest periodictable future https://github.com/bjodah/pycompilation/archive/v0.2.21.tar.gz https://github.com/sympy/sympy/archive/master.zip
python -c "import sympy; print(sympy.__version__)"
python -c "import pycompilation; print(pycompilation.__version__)"
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
python setup.py build_ext -i
PYTHONPATH=.:$PYTHONPATH py.test --slow
