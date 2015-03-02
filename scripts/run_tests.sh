#!/bin/bash
# Usage, e.g.:
# ./scripts/run_tests.sh --ignore examples/

export PKG_NAME=chemreac

PYTEST_ARGS=()

# py.test might use either 'python' or 'python3'
PYTHON_EXE=$(head -1 $(which py.test) | cut -f2 -d!)

# Check dependencies for the py.test test command below,
# note that the package itself might depend on more packages
# (see ../requirements.txt)
if ! $PYTHON_EXE -c "import pytest" > /dev/null 2>&1; then
    >&2 echo "Error, could not import pytest, please install pytest."
fi

if ! $PYTHON_EXE -c "import pytest_pep8" > /dev/null 2>&1; then
    echo "Could not import pytest_pep8, install pytest-pep8 if you want it."
else
    PYTEST_ARGS+=(--pep8)
fi

if ! $PYTHON_EXE -c "import pytest_cov" > /dev/null 2>&1; then
    echo "Could not import pytest_cov, install pytest-cov if you want it."
else
    PYTEST_ARGS+=(--cov $PKG_NAME --cov-report html)
fi

if ! $PYTHON_EXE -c "import $PKG_NAME" > /dev/null 2>&1; then
    if ! $PYTHON_EXE -c "import pycompilation" > /dev/null 2>&1; then
        >&2 echo "Error, could not import pycompilation, please install pycompilation."
    else
        >&2 echo "$PKG_NAME has not been built, building inplace..."
        $PYTHON_EXE setup.py build_ext --inplace
        export TEST_INSTALL=1
    fi
fi
if [[ "$TEST_INSTALL" != "1" ]]; then
    # Extract absolute path of script, from:
    # http://unix.stackexchange.com/a/9546
    absolute_repo_path_x="$(readlink -fn -- "$(dirname $0)/.."; echo x)"
    absolute_repo_path="${absolute_repo_path_x%x}"
    export PYTHONPATH=$absolute_repo_path:$PYTHONPATH
fi
echo "About to run the full test suite. It can take several minutes..."
set -xe  # bash: echo commands, exit on failure
py.test --slow --veryslow --doctest-modules ${PYTEST_ARGS[@]} $@
