#!/bin/bash -e
PY_VERSION=$1
MINICONDA_PATH=$2
CONDA_VERSION=$3

if [[ $(uname -m) == "x86_64" ]]; then
    # http://stackoverflow.com/a/106399
    ARCH="x86_64"
else
    ARCH="x86"
fi

if [[ "${PY_VERSION:0:1}" == "2" ]]; then
    if [[ "${PY_VERSION}" != "2.7" ]]; then
        echo "PY_VERSION was not set to 2.7"
        exit 1
    fi
    wget --quiet http://repo.continuum.io/miniconda/Miniconda-${CONDA_VERSION}-Linux-${ARCH}.sh -O miniconda.sh;
else
    if [[ "${PY_VERSION:0:1}" == "3.4" ]]; then
        echo "PY_VERSION was not set to 2.7 or 3.4"
        exit 1
    fi
    wget --quiet http://repo.continuum.io/miniconda/Miniconda3-${CONDA_VERSION}-Linux-${ARCH}.sh -O miniconda.sh;
fi
bash miniconda.sh -b -p $MINICONDA_PATH
rm miniconda.sh
export PATH="$MINICONDA_PATH/bin:$PATH"
hash -r
