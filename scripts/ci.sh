#!/bin/bash -x
PKG_NAME=${1:-${CI_REPO##*/}}
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${1^^}_RELEASE_VERSION=\$CI_BRANCH
    echo ${CI_BRANCH} | tail -c +2 > __conda_version__.txt
fi

set -e
for PY in python2 python3; do
    $PY -m pip install --user -e .[all]
    if [[ $PY == *2 ]]; then
        PYTHON=$PY ./scripts/run_tests.sh
    else
        PYTHON=$PY ./scripts/run_tests.sh ${@:2}
    fi
done

[[ $(python setup.py --version) =~ ^[0-9]+.[0-9]+* ]]  # make sure ./setup.py --version returns X.X*

python2 setup.py sdist
cp dist/${PKG_NAME}-*.tar.gz /tmp
(cd /; python2 -m pip install --force-reinstall /tmp/${PKG_NAME}-*.tar.gz; python2 -c "import $PKG_NAME")

# Make sure repo is pip installable from git-archive zip
git archive -o /tmp/$PKG_NAME.zip HEAD
(cd /; python3 -m pip install --force-reinstall /tmp/$PKG_NAME.zip; python3 -c "import ${PKG_NAME}")
