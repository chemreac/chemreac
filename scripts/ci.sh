#!/bin/bash -xu
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${1^^}_RELEASE_VERSION=\$CI_BRANCH
    echo ${CI_BRANCH} | tail -c +2 > __conda_version__.txt
fi
rm -r build/ dist/* */*.so
set -e
python2.7 setup.py sdist
python2.7 -m pip install --user -e .[all]
python3.4 -m pip install --user -e .[all]
PYTHONPATH=$(pwd) PYTHON=python2.7 ./scripts/run_tests.sh
PYTHONPATH=$(pwd) PYTHON=python3.4 ./scripts/run_tests.sh ${@:2}
! grep "DO-NOT-MERGE!" -R . --exclude ci.sh
