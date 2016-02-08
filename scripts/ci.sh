#!/bin/bash -xu
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${1^^}_RELEASE_VERSION=\$CI_BRANCH
    echo ${CI_BRANCH} | tail -c +2 > __conda_version__.txt
fi
rm -r build/ dist/* */*.so
set -e
python2 setup.py sdist
pip install dist/*.tar.gz
(cd /; python2.7 -m pytest --pyargs $1)
pip3 install dist/*.tar.gz
(cd /; python3 -m pytest --pyargs $1)
pip install --user -e .[all]
pip3 install --user -e .[all]
PYTHONPATH=$(pwd) ./scripts/run_tests.sh $@
! grep "DO-NOT-MERGE!" -R . --exclude ci.sh
