#!/bin/bash -xu
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${1^^}_RELEASE_VERSION=\$CI_BRANCH
    echo ${CI_BRANCH} | tail -c +2 > __conda_version__.txt
fi
rm -r build/ dist/* */*.so
set +e
python2 setup.py sdist
pip install dist/*.tar.gz
(cd /; python2.7 -m pytest --pyargs $1)
pip3 install dist/*.tar.gz
(cd /; python3 -m pytest --pyargs $1)
python2 setup.py build_ext -i
python3 setup.py build_ext -i
PYTHONPATH=$(pwd) ./scripts/run_tests.sh $@
./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg
! grep "DO-NOT-MERGE!" -R . --exclude ci.sh
