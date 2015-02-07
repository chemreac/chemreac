#!/bin/bash -x
if [[ "$TRAVIS_REPO_SLUG" == "${GITHUB_USER}/${GITHUB_REPO}" ]] && [[ "$TRAVIS_PULL_REQUEST" == "false" ]]; then
    if [[ "$TRAVIS_BRANCH" == "master" ]]; then
        if [[ "$BUILD_DOCS" == "1" ]]; then
            # Build the documentation
            echo -e "Building docs...\n"
            set -x # Verbose
            ./scripts/build_docs.sh

            echo -e "Publishing pages...\n"
            WORKDIR=`pwd`
            cd $HOME
            git config --global user.email "travis@travis-ci.org"
            git config --global user.name "travis-ci"
            set +x # Silent (protect GH_TOKEN)
            echo "Cloning github repo: ${TRAVIS_REPO_SLUG}"
            git clone --quiet https://${GH_TOKEN}@github.com/${GITHUB_USER}/${GITHUB_USER}.github.io ${GITHUB_USER}.github.io > /dev/null
            set -x # Verbose
            cd ${GITHUB_USER}.github.io/
            git branch -D master
            git checkout --orphan master
            git rm -rf . > /dev/null
            cp -R ${WORKDIR}/gh-pages-skeleton/* .
            cp ${WORKDIR}/gh-pages-skeleton/.* .
            cp -R ${WORKDIR}/docs/_build/html ./docs
            cp -R ${WORKDIR}/htmlcov .
            git add -f . > /dev/null
            git commit -m "Lastest pages from successful travis build $TRAVIS_BUILD_NUMBER"
            set +x # Silent (protect GH_TOKEN)
            git push -f origin master >/dev/null 2>&1
            set -x
            echo -e "...published to ${GITHUB_USER}.github.io\n"
        fi
    fi
    if [[ "$CHEMREAC_RELEASE_VERSION" == v* ]]; then  # version-tagged release => upload to binstar
        conda install binstar
        cat __conda_version__.txt # DEBUGGING failed upload
        ls /home/travis/miniconda/conda-bld/linux-64/ # DEBUGGING failed upload
        export MY_CONDA_PKG=$(conda build --output conda-recipe | tail -n 1)
        set +x  # Silent (protect token in Travis log)
        binstar -t $BINSTAR_TOKEN upload --force ${MY_CONDA_PKG/--/-${CHEMREAC_RELEASE_VERSION#v}-}
        set -x
    fi
fi
