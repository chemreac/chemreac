#!/bin/bash -e
# Usage
#   $ ./scripts/run_tests.sh
# or
#   $ ./scripts/run_tests.sh --cov pycvodes --cov-report html
export MPLBACKEND=Agg
ulimit -v 2048000
${PYTHON:-python3} -m pytest --doctest-modules --pep8 --flakes --slow --veryslow $@
${PYTHON:-python3} -m doctest README.rst
