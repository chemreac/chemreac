#!/bin/bash
#!/bin/bash
# Usage, e.g.:
# ./scripts/run_tests.sh --ignore examples/

export PKG_NAME=chemreac

# Check dependencies for the py.test test command below,
# note that the package itself might depend on more packages
# (see ../requirements.txt)
if ! python -c "import pytest" > /dev/null 2>&1; then
    >&2 echo "Error, could not import pytest, please install pytest."
fi
if ! python -c "import pytest_pep8" > /dev/null 2>&1; then
    >&2 echo "Error, could not import pytest_pep8, please install pytest-pep8."
fi
if ! python -c "import pytest_cov" > /dev/null 2>&1; then
    >&2 echo "Error, could not import pytest_cov, please install pytest-cov."
fi
if ! python -c "import $PKG_NAME" > /dev/null 2>&1; then
    if ! python -c "import pycompilation" > /dev/null 2>&1; then
        >&2 echo "Error, could not import pycompilation, please install pycompilation."
    else
        >&2 echo "$PKG_NAME has not been built, building inplace..."
        python setup.py build_ext --inplace
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
set -xe  # bash: echo commands, exit on failure
py.test --slow --veryslow --pep8 --doctest-modules --cov $PKG_NAME --cov-report html $@
