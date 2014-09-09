#!/bin/bash

PY_VERSION=$1
MINICONDA_PATH=$2

ENV_NAME=test-env
export PATH=$MINICONDA_PATH/bin:$PATH

conda create --quiet -n $ENV_NAME python=${PY_VERSION} numpy=1.8.1 scipy=0.14 matplotlib=1.4.0 cython=0.20.2 pip
source activate $ENV_NAME
pip install --quiet argh mako quantities pytest pytest-pep8 sphinx numpydoc sphinx_rtd_theme periodictable future https://github.com/bjodah/pycompilation/archive/v0.3.3.tar.gz https://github.com/bjodah/pycodeexport/archive/v0.0.4.tar.gz https://github.com/sympy/sympy/archive/master.zip
