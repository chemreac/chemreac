#!/bin/bash -xeu
# Usage:
#
#    $ ./scripts/release.sh v1.2.3 ~/anaconda2/bin myserver.example.com GITHUB_USER GITHUB_REPO upstream
#

if [[ $1 != v* ]]; then
    echo "Argument does not start with 'v'"
    exit 1
fi
VERSION=${1#v}
CONDA_PATH=$2
SERVER=$3
GITHUB_USER=$4
GITHUB_REPO=$5
REMOTE=$6
find . -type f -iname "*.pyc" -exec rm {} +
find . -type f -iname "*.o" -exec rm {} +
find . -type f -iname "*.so" -exec rm {} +
find . -type d -name "__pycache__" -exec rmdir {} +
find . -type f -iname ".coverage.*" -exec rm {} +
./scripts/check_clean_repo_on_master.sh
for DIR in build/ dist/ docs/_build/ *.egg-info/ .cache/; do
    if [[ -e $DIR ]]; then
        rm -r $DIR
    fi
done
./scripts/check_clean_repo_on_master.sh
cd $(dirname $0)/..
# PKG will be name of the directory one level up containing "__init__.py" 
PKG=$(find . -maxdepth 2 -name __init__.py -print0 | xargs -0 -n1 dirname | xargs basename)
PKG_UPPER=$(echo $PKG | tr '[:lower:]' '[:upper:]')
./scripts/run_tests.sh
env ${PKG_UPPER}_RELEASE_VERSION=v$VERSION python setup.py sdist
if [[ -e ./scripts/generate_docs.sh ]]; then
    env ${PKG_UPPER}_RELEASE_VERSION=v$VERSION ./scripts/generate_docs.sh
fi
for CONDA_PY in 2.7 3.4 3.5; do
    continue  # we build the conda recipe on another host for now..
    for CONDA_NPY in 1.11; do
        PATH=$2:$PATH ./scripts/build_conda_recipe.sh v$VERSION --python $CONDA_PY --numpy $CONDA_NPY
    done
done

# All went well, add a tag and push it.
git tag -a $1 -m $1
git push $REMOTE
git push $REMOTE --tags
twine upload dist/${PKG}-$VERSION.tar.gz

set +x
echo ""
echo "    You may now create a new github release at with the tag \"v$VERSION\" and name"
echo "    it \"${PKG}-${VERSION}\". Here is a link:"
echo "        https://github.com/${GITHUB_USER}/${GITHUB_REPO}/releases/new "
echo "    name the release \"${PKG}-${VERSION}\", and don't foreget to manually attach the file:"
echo "        $(openssl sha256 $(pwd)/dist/${PKG}-${VERSION}.tar.gz)"
echo "    Then run:"
echo ""
echo "        $ ./scripts/post_release.sh $1 $SERVER $4"
echo ""
