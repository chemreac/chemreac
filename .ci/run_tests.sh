#!/bin/bash -x
PKG_NAME=${1:-${CI_REPO##*/}}
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${1^^}_RELEASE_VERSION=\$CI_BRANCH
fi

python3 -m pip install --user -e .[all]

set -e

(cd tests-native; make -B CONTEXT=valgrind EXTRA_COMPILE_ARGS='-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC' test)
(cd tests-native; make -B CXX=clang++-8 CC=clang-8 OPTIMIZE=1 WITH_OPENMP=0 EXTRA_COMPILE_ARGS='-fsanitize=address -DNDEBUG' test)

rm -r build/
CFLAGS="-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" python3 setup.py build_ext -i
bash -c "ulimit -v 2048000; ./scripts/run_tests.sh"

rm -r build/
CC=clang-8 \
  CXX=clang++-8 \
  CFLAGS="-fsanitize=address -UNDEBUG" \
  python3 setup.py build_ext -i
ASAN_OPTIONS=detect_leaks=0 LD_PRELOAD=/usr/lib/llvm-8/lib/clang/8.0.1/lib/linux/libclang_rt.asan-x86_64.so ./scripts/run_tests.sh "${@:2}"

rm -r build/
python3 setup.py sdist
cp dist/${PKG_NAME}-*.tar.gz /tmp
(cd /; python3 -m pip uninstall -y chemreac; python3 -m pip install /tmp/${PKG_NAME}-*.tar.gz; python3 -c "import $PKG_NAME")

# Make sure repo is pip installable from git-archive zip
git archive -o /tmp/$PKG_NAME.zip HEAD
(cd /; python3 -m pip uninstall -y chemreac; python3 -m pip install /tmp/$PKG_NAME.zip; python3 -c "import ${PKG_NAME}")
