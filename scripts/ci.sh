#!/bin/bash -x
PKG_NAME=${1:-${CI_REPO##*/}}
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${1^^}_RELEASE_VERSION=\$CI_BRANCH
    echo ${CI_BRANCH} | tail -c +2 > __conda_version__.txt
fi

set -e
python3 -m pip install --user -e .[all]
./scripts/run_tests.sh ${@:2}
python3 -m pip uninstall -y $PKG_NAME

python3 setup.py sdist
cp dist/${PKG_NAME}-*.tar.gz /tmp
(cd /; python3 -m pip install --force-reinstall /tmp/${PKG_NAME}-*.tar.gz; python3 -c "import $PKG_NAME")
python3 -m pip uninstall -y $PKG_NAME

# Make sure repo is pip installable from git-archive zip
git archive -o /tmp/$PKG_NAME.zip HEAD
(cd /; python3 -m pip install --force-reinstall /tmp/$PKG_NAME.zip; python3 -c "import ${PKG_NAME}")
