#!/bin/bash

PY_VERSION=$1
ENV_NAME=$2

conda create --quiet -n $ENV_NAME python=${PY_VERSION} numpy=1.9.2 scipy=0.16.0 matplotlib=1.4.3 cython=0.22.1 mako periodictable quantities pytest future pip sphinx numpydoc pycompilation pycodeexport sympy
source activate $ENV_NAME

pip install --quiet argh pytest-pep8 pytest-cov python-coveralls sphinx_rtd_theme mpld3
# mako quantities periodictable
