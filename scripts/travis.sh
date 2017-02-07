#!/bin/bash
conda config --add channels conda-forge
conda config --add channels bjodah
conda config --set always_yes yes
conda install conda-build
cp -r conda-recipe/ $1
cd /tmp/conda-recipe
sed -i -E -e "s~git_url:(.+)~path: $TRAVIS_BUILD_DIR~" $1/meta.yaml
