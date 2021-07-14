#!/bin/bash -x
PKG_NAME=${1:-${CI_REPO##*/}}
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${1^^}_RELEASE_VERSION=\$CI_BRANCH
fi

python3 -m pip install --user -e .[all]
python3 setup.py build_ext -i
mkdir -p dist
cp -r chemreac dist/.
ls dist/chemreac
rm -r build
set -e
( cd dist; python3 -c "import chemreac" )

(cd tests-native; make -B CONTEXT=valgrind EXTRA_COMPILE_ARGS='-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC' test)
(cd tests-native; make -B CXX=clang++-12 CC=clang-11 OPTIMIZE=1 WITH_OPENMP=0 EXTRA_COMPILE_ARGS='-fsanitize=address -DNDEBUG' test)

CFLAGS="-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC $CFLAGS" python3 setup.py build_ext -i
gdb -q -ex "set confirm off" -ex r -ex bt -ex q -args pythno3 -m pytest
exit 1  # DO-NOT-MERGE (gdb debugging via CI server)
#bash -c "ulimit -v 2048000; ./scripts/run_tests.sh"

rm -rf build/
CC=clang-11 \
  CXX=clang++-12 \
  CFLAGS="-fsanitize=address -UNDEBUG $CFLAGS" \
  python3 setup.py build_ext -i
PYTHONMALLOC=malloc ASAN_OPTIONS=detect_leaks=0 LD_PRELOAD=$(clang++-12 -print-file-name libclang_rt.asan-$(uname -m).so) ./scripts/run_tests.sh "${@:2}"

python3 -m pip uninstall -y ${PKG_NAME}

rm -rf build/

python3 setup.py sdist
cp dist/${PKG_NAME}-*.tar.gz /tmp
(mkdir /tmp/sdist_tar_gz; cd /tmp/sdist_tar_gz; tar xf ../${PKG_NAME}-*.tar.gz; cd ${PKG_NAME}-*/; python3 setup.py build_ext -i; PYTHONPATH=$(pwd); python3 -c "import $PKG_NAME")

# Make sure repo is pip installable from git-archive zip
(mkdir /tmp/archive_zip; cd /tmp/archive_zip; unzip -q ../HEAD.zip; python3 setup.py build_ext -i; PYTHONPATH=$(pwd); python3 -c "import $PKG_NAME")
