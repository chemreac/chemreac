#!/bin/bash

PYTHON_VERSION=$1
ENV_NAME=$2

conda create --quiet -n $ENV_NAME python=${PYTHON_VERSION} numpy=1.9.2 scipy=0.16.0 matplotlib=1.4.3 cython=0.22.1 mako periodictable quantities pytest future sphinx numpydoc pycompilation pycodeexport sympy sphinx_rtd_theme argh pytest-pep8 pytest-cov python-coveralls mpld3
