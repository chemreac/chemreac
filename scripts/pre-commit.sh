#!/bin/bash
if git diff-index --quiet HEAD --; then
    # no changes between index and working copy; just run tests
    OMP_NUM_THREADS=1 python -m pytest --pep8 --maxfail 1
    RESULT=$?
else
    # Test the version that's about to be committed,
    # stashing all unindexed changes
    git stash -q --keep-index
    py.test --pep8 --maxfail 1
    RESULT=$?
    git stash pop -q
fi

if [[ "$RESULT" -ne "0" ]]; then
    echo "Aborting commit.  Fix above errors or do 'git commit --no-verify'."
    exit 1
fi
