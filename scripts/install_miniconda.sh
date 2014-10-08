#!/bin/bash -e
PY_VERSION=$1
INSTALL_PATH=$2

if [[ "${PY_VERSION:0:1}" == "2" ]]; then
    if [[ "${PY_VERSION:0:1}" == "2.7" ]]; then
        echo "PY_VERSION was not set to 2.7"
        exit 1
    fi
    # This saves us some downloading for this version
    wget --quiet http://repo.continuum.io/miniconda/Miniconda-3.7.0-Linux-x86_64.sh -O miniconda.sh;
else
    if [[ "${PY_VERSION:0:1}" == "3.4" ]]; then
        echo "PY_VERSION was not set to 2.7 or 3.4"
        exit 1
    fi
    wget --quiet http://repo.continuum.io/miniconda/Miniconda3-3.7.0-Linux-x86_64.sh -O miniconda.sh;
fi
if [[ $? != 0 ]]; then
    # Only matters when -e not used
    echo "Failed to get Miniconda."
    exit 1
fi

chmod +x miniconda.sh
./miniconda.sh -b -p $INSTALL_PATH
export PATH="$INSTALL_PATH/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update --quiet conda
conda info -a
