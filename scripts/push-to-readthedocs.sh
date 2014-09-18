#!/bin/bash -x
GH_USER=bjodah
GH_REPO=chemreac
if [ "$TRAVIS_PYTHON_VERSION" == "3.4" ] && [ "$TRAVIS_REPO_SLUG" == "${GH_USER}/${GH_REPO}" ] && [ "$TRAVIS_PULL_REQUEST" == "false" ] && [ "$TRAVIS_BRANCH" == "master" ]; then
    echo -e "Triggering readthedocs webhook...\n"
    curl -X POST http://readthedocs.org/build/chemreac
fi
