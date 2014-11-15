#!/bin/bash -x
if [ "$TRAVIS_REPO_SLUG" == "${GITHUB_USER}/${GITHUB_REPO}" ] && [ "$TRAVIS_PULL_REQUEST" == "false" ] && [ "$TRAVIS_BRANCH" == "master" ]; then
    # Publish binary if version does not end with .dev
    export DISTUTILS_DEBUG=  # less verbose setup.py --version
    export VERSION=$(python setup.py --version)
    if [[ "${VERSION}" != *.dev ]]; then
        conda install binstar
        export MY_CONDA_PKG=$(conda build --output conda-recipe | tail -n 1)
        set +x # Silent (protect token in Travis log)
        binstar -t $BINSTAR_TOKEN upload --force ${MY_CONDA_PKG/--/-$VERSION-}
        set -x
    fi

    if [ "$BUILD_DOCS" == "1" ]; then
        # Build the documentation
        echo -e "Building docs...\n"
        set -x # Verbose
        ./scripts/build_docs.sh

        echo -e "Publishing docs...\n"
        WORKDIR=`pwd`
        cd $HOME
        git config --global user.email "travis@travis-ci.org"
        git config --global user.name "travis-ci"
        set +x # Silent (protect GH_TOKEN)
        echo "Cloning github repo: ${TRAVIS_REPO_SLUG}"
        git clone --quiet https://${GH_TOKEN}@github.com/${TRAVIS_REPO_SLUG} ${GITHUB_REPO} > /dev/null
        set -x # Verbose
        cd ${GITHUB_REPO}
        git branch -D gh-pages
        git checkout --orphan gh-pages
        git rm -rf . > /dev/null
        cp -R ${WORKDIR}/gh-pages-skeleton/* .
        cp ${WORKDIR}/gh-pages-skeleton/.* .
        cp -R ${WORKDIR}/docs/_build/html ./docs
        cp -R ${WORKDIR}/htmlcov .
        git add -f . > /dev/null
        git commit -m "Lastest docs from successful travis build $TRAVIS_BUILD_NUMBER"
        git push -f origin gh-pages
        echo -e "Published docs to gh-pages.\n"
    fi
fi
