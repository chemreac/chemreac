#!/bin/bash -e
# this script assumes conda is in $PATH

export PYTHON_VERSION=$1
export CONDA_PY=$2
export ENV_NAME=$3
export RUN_TESTS=$4
export SETUP_TEST_ENV=${5:-1}

export CONDA_PATH=$(conda info --system | grep sys.prefix | cut -d: -f2 | sed -e 's/^ *//')
if [[ "$SETUP_TEST_ENV" == "1" ]]; then
    source ./scripts/setup_conda_testenv.sh $PYTHON_VERSION $ENV_NAME
    export LIBRARY_PATH=$CONDA_PATH/envs/$ENV_NAME/lib:$LIBRARY_PATH
else
    export LIBRARY_PATH=$CONDA_PATH/lib:$LIBRARY_PATH
fi
conda info
python --version
python setup.py --version
export DISTUTILS_DEBUG=1
conda build --quiet conda-recipe
conda install --quiet chemreac --use-local
if [[ "$RUN_TESTS" == "1" ]]; then
    ./scripts/run_tests.sh
fi
