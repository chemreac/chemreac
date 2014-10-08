#!/bin/bash -x
GH_USER=bjodah
GH_REPO=chemreac
if [ "$TRAVIS_PYTHON_VERSION" == "3.4" ] && [ "$TRAVIS_REPO_SLUG" == "${GH_USER}/${GH_REPO}" ] && [ "$TRAVIS_PULL_REQUEST" == "false" ] && [ "$TRAVIS_BRANCH" == "master" ]; then
    # ReadTheDocs
    echo -e "Triggering readthedocs webhook...\n"
    curl -X POST http://readthedocs.org/build/chemreac

    # Github pages
    echo -e "Building docs...\n"
    ./scripts/build_docs.sh

    echo -e "Publishing docs...\n"
    cd $HOME
    git config --global user.email "travis@travis-ci.org"
    git config --global user.name "travis-ci"
    git clone --quiet https://git@github.com/${GH_USER}/${GH_REPO} ${GH_REPO} > /dev/null
    cd ${GH_REPO}
    git branch -D gh-pages
    git checkout --orphan gh-pages
    git rm -rf .
    cp -R ${WORKDIR}/gh-pages-skeleton/* .
    cp ${WORKDIR}/gh-pages-skeleton/.* .
    cp -R ${WORKDIR}/docs/_build/html ./docs
    cp -R ${WORKDIR}/htmlcov .
    git add -f .
    git commit -m "Lastest docs on successful travis build $TRAVIS_BUILD_NUMBER auto-pushed to gh-pages"
    git push -f origin gh-pages
    echo -e "Published docs to gh-pages.\n"

    # Trigger drone.io build (more extensive testing)
    curl -X POST "https://drone.io/hook?id=github.com/bjodah/chemreac&token=LSWwdhBX47ZnYPs4gtrX"

    # Push wheel to binstar.org
    conda config --add channels http://conda.binstar.org/USER
    conda build conda-recipe/
    binstar -t `cat $HOME/binstar_token` upload --force *.tar.bz2 # bash -x workaround
fi
