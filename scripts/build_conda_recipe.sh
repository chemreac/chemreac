#!/bin/bash -ex
# Usage:
#
#    $ ./scripts/build_conda_recipe.sh v1.2.3 --python 2.7
#
echo ${1#v}>__conda_version__.txt
cleanup() {
    rm __conda_version__.txt
}
trap cleanup INT TERM EXIT

if [[ -d build/ ]]; then
    rm -r build/
fi
conda build ${@:2} ./conda-recipe/
