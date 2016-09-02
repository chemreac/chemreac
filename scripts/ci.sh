#!/bin/bash -xu
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${1^^}_RELEASE_VERSION=\$CI_BRANCH
    echo ${CI_BRANCH} | tail -c +2 > __conda_version__.txt
fi
rm -r build/ dist/* */*.so
set -e
python2.7 setup.py sdist
for PY in python2 python3; do
    $PY -m pip install --user -e .[all]
    if [[ $PY == *2 ]]; then
        PYTHON=$PY ./scripts/run_tests.sh
    else
        PYTHON=$PY ./scripts/run_tests.sh ${@:2}
    fi
done
! grep "DO-NOT-MERGE!" -R . --exclude ci.sh
